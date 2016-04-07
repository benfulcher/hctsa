function [ifeat, testStat, testStat_rand] = TS_TopFeatures(whatData,whatTestStat,doNull,varargin)
% TS_TopFeatures    Top individual features for discriminating labeled time series
%
% This function compares each feature in an hctsa dataset individually for its
% ability to separate the labeled classes of time series according to a given
% test statistic.
%
% Can also compare this performance to a set of randomized null features to
% evaluate the statistical significance of the result.
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData, default: 'raw')
% whatTestStat, the test statistic to quantify the goodness of each feature
%               (one of: 'tstat', 'svm', 'linear', 'diaglinear')
% doNull, (binary) whether to compute a null distribution using permutations of the class
%           labels
%
% ***Additional plotting options***:
% 'whatPlots', can specify what output plots to produce (cell of strings), e.g.,
%               specify {'histogram','distributions','cluster'} to produce all
%               three possible output plots (this is the default).
% 'numTopFeatures', can specify the number of top features to analyze, both in
%                   terms of the list of outputs, the histogram plots, and the
%                   cluster plot.
% 'numHistogramFeatures', can optionally also set a custom number of histograms
%                       to display (often want to set this lower to avoid producing
%                       large numbers of figures)
%
%---EXAMPLE USAGE:
%
% TS_TopFeatures('norm','tstat',1,'whatPlots',{'histogram','distributions','cluster'},'numTopFeatures',20);
%
%---OUTPUTS:
% ifeat, the ordering of operations by their performance
% testStat, the test statistic (whatTestStat) for each operation
% testStat_rand, test statistics making up the null distributions

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% --------------------------------------------------------------------------
%% Check inputs
% --------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'raw';
end
if nargin < 2 || isempty(whatTestStat)
    whatTestStat = 'fast_linear'; % Way faster than proper prediction models
    fprintf(1,'Using ''%s'' test statistic by default\n', whatTestStat);
end
if nargin < 3
    doNull = 1; % compute an empirical null distribution by randomizing class labels
end

% Use an inputParser to control plotting options:
inputP = inputParser;
default_whatPlots = {'histogram','distributions','cluster'};
check_whatPlots = @(x) iscell(x) || ischar(x);
addParameter(inputP,'whatPlots',default_whatPlots,check_whatPlots);
default_numTopFeatures = 40;
addParameter(inputP,'numTopFeatures',default_numTopFeatures,@isnumeric);
default_numHistogramFeatures = 10;
addParameter(inputP,'numHistogramFeatures',default_numHistogramFeatures,@isnumeric);
parse(inputP,varargin{:});

whatPlots = inputP.Results.whatPlots;
numTopFeatures = inputP.Results.numTopFeatures;
numHistogramFeatures = inputP.Results.numHistogramFeatures;
clear inputP;

% --------------------------------------------------------------------------
%% Load the data
% --------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
numOps = length(Operations);

%-------------------------------------------------------------------------------
%% Check that grouping information exists:
%-------------------------------------------------------------------------------
if ~isfield(TimeSeries,'Group')
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
timeSeriesGroup = [TimeSeries.Group]'; % Use group form
numGroups = length(unique(timeSeriesGroup));
groupNames = load(whatDataFile,'groupNames');
groupNames = groupNames.groupNames;

% --------------------------------------------------------------------------
%% Define the train/test classification rate function, fn_testStat
% --------------------------------------------------------------------------
% Also the chanceLine -- where you'd expect by chance (for equiprobable groups...)
if ismember(whatTestStat,{'linear','linclass','fast_linear','diaglinear','svm','svm_linear'})
    % Percentage correct:
    fn_testStat = GiveMeFunctionHandle(whatTestStat,numGroups);
    % fn_testStat = @(XTrain,yTrain,XTest,yTest) ...
    %         GiveMeCfn(whatTestStat,XTrain,yTrain,XTest,yTest,numGroups,0)*100;
    chanceLine = 100/numGroups;
    cfnUnit = '%';
end

switch whatTestStat
case {'linear','linclass','fast_linear'}
    fprintf(1,'Using a linear classifier\n');
    cfnName = 'linear classifier';
    % cfnUnit = '%';
    % fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
    %                 mean(yTest == classify(Xtest,XTrain,yTrain,'linear'))*100;
    % chanceLine = 100/length(unique(timeSeriesGroup));
case 'diaglinear'
    fprintf(1,'A Naive Bayes classifier\n');
    cfnName = 'naive bayes classifier';
    % cfnUnit = '%';
    % fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
    %                 mean(yTest == classify(Xtest,XTrain,yTrain,'diaglinear'))*100;
    % chanceLine = 100/length(unique(timeSeriesGroup));
case {'svm','svmlinear'}
    fprintf(1,'A linear support vector machine\n');
    cfnName = 'SVM classifier';
    % cfnUnit = '%';
    % fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
    %                 mean(yTest == predict(fitcsvm(XTrain,yTrain, ...
    %                         'KernelFunction','linear'),Xtest))*100;
    % chanceLine = 100/length(unique(timeSeriesGroup));
    if numGroups > 2
        error('SVM not supported for multiclass classification with more than 2 classes');
    end
case {'ttest','tstat'}
    fprintf(1,'A Welch''s t-statistic\n');
    cfnName = 'Welch t-stat';
    cfnUnit = '';
    if numGroups > 2
        error('Cannot use t-test as test statistic with more than two groups :/');
    end
    fn_testStat = @(XTrain,yTrain,Xtest,yTest) fn_tStat(XTrain(yTrain==1),XTrain(yTrain==2));
    chanceLine = 0; % by chance, t-stats averge to zero
otherwise
    error('Unknown method ''%s''',whatTestStat)
end


% --------------------------------------------------------------------------
%%                     Loop over all features
% --------------------------------------------------------------------------
% Use the same data for training and testing:
fprintf(1,'Comparing the (in-sample) performance of %u operations...',length(Operations));
timer = tic;
testStat = giveMeStats(TS_DataMat,timeSeriesGroup);
fprintf(1,' Done in %s.\n',BF_thetime(toc(timer)));

% Give mean and that expected from random classifier (may be a little overfitting)
fprintf(1,['Mean %s across %u operations = %4.2f%s\n' ...
            '(Random guessing for %u equiprobable classes = %4.2f%s)\n'], ...
        cfnName,numOps,nanmean(testStat),cfnUnit,numGroups,chanceLine,cfnUnit);

% --------------------------------------------------------------------------
%% Display information about the top topN operations
% --------------------------------------------------------------------------
[testStat_sort, ifeat] = sort(testStat,'descend'); % bigger is better

isNaN = isnan(testStat_sort);
testStat_sort = testStat_sort(~isNaN);
ifeat = ifeat(~isNaN);

% List the top features:
topN = min(numHistogramFeatures,length(Operations));
for i = 1:topN
    fprintf(1,'[%u] %s (%s) -- %4.2f%%\n',Operations(ifeat(i)).ID, ...
            Operations(ifeat(i)).Name,Operations(ifeat(i)).Keywords,testStat_sort(i));
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%% Plot outputs
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Histogram of distribution of test statistics for labeled and null data
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'histogram'))
    % A figure to show the distribution of test statistics across all
    % features:

    %-------------------------------------------------------------------------------
    %% Compute null distribution
    %-------------------------------------------------------------------------------
    numRepeats = 10;
    testStat_rand = zeros(numOps,numRepeats);
    if doNull
        fprintf(1,'Now for %u nulls... ',numRepeats);
        tic
        for j = 1:numRepeats
            if j<numRepeats
                fprintf(1,'%u,',j);
            else
                fprintf(1,'%u',j);
            end
            % Shuffle labels:
            groupLabels = timeSeriesGroup(randperm(length(timeSeriesGroup)));
            testStat_rand(:,j) = giveMeStats(TS_DataMat,groupLabels);
        end
        fprintf(1,'\n%u %s statistics computed in %s.\n',numOps*numRepeats,...
                                        cfnName,BF_thetime(toc(timer)));
    end

    f = figure('color','w'); hold on
    colors = BF_getcmap('spectral',5,1);
    if ~doNull
        h_real = histogram(testStat,'Normalization','probability',...
                    'BinMethod','auto','EdgeColor',colors{5},'FaceAlpha',0);
        maxH = max(h_real.Values);
    else
        % Plot both real distribution, and null distribution:
        numBins = 20;
        allTestStat = [testStat(:);testStat_rand(:)];
        minMax = [min(allTestStat),max(allTestStat)];
        histEdges = linspace(minMax(1),minMax(2),numBins+1);
        h_null = histogram(testStat_rand(:),histEdges,'Normalization','probability','FaceColor',colors{1});
        h_real = histogram(testStat,histEdges,'Normalization','probability','EdgeColor',colors{5});
        maxH = max([max(h_real.Values),max(h_null.Values)]);
        plot(mean(testStat_rand(:))*ones(2,1),[0,maxH],'--','color',colors{1},'LineWidth',2)
    end

    % Add chance line:
    l_chance = plot(chanceLine*ones(2,1),[0,maxH],'--k');
    % Add mean of real distribution:
    l_mean = plot(mean(testStat)*ones(2,1),[0,maxH],'--','color',colors{5},'LineWidth',2);

    % Labels:
    xlabel(sprintf('Individual %s across %u features',cfnName,numOps))
    ylabel('Probability')

    % Legend:
    if doNull,
        legend([h_null,h_real,l_chance,l_mean],'null','real','chance','real mean')
    else
        legend([h_real,l_chance,l_mean],{'real','chance','real mean'});
    end
end

%-------------------------------------------------------------------------------
% Distributions across classes for top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'distributions'))
    subPerFig = 5; % subplots per figure
    ks_or_hist = 'ks'; % view as either histograms or kernel-smoothed distributions

    % Set the colors to be assigned to groups:
    if numGroups<=5
        colors = BF_getcmap('set1',5,1);
        if numGroups==2, colors = colors([2,4]); end
    else
        colors = BF_getcmap('dark2',numGroups,1);
        if length(colors) < numGroups
            % Too many groups for a custom colormap, just space them along the jet colormap:
            colors = jet(numGroups);
            colors = arrayfun(@(x)colors(x,:),1:size(colors,1),'UniformOutput',0);
        end
    end

    % Space the figures out properly:
    numFigs = ceil(numTopFeatures/subPerFig);

    for figi = 1:numFigs
        if figi*subPerFig > length(ifeat)
            break % We've exceeded number of features
        end
        r = ((figi-1)*subPerFig+1:figi*subPerFig);
        featHere = ifeat(r); % features to plot on this figure

        f = figure('color','w');
        f.Position(3:4) = [588, 612]; % make longer
        for opi = 1:subPerFig
            op_ind = featHere(opi);
            subplot(subPerFig,1,opi); box('on'); hold on

            switch ks_or_hist
            case 'ks'
                % Plot distributions first for the sake of the legend
                linePlots = cell(numGroups,1);
                for i = 1:numGroups
                    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                    [f,x,linePlots{i}] = BF_plot_ks(featVector,colors{i},0,2,20,1);
                    % [f, x] = ksdensity(featVector);
                    % % Plot only the range represented in the actual feature vector
                    % rMatch = (arrayfun(@(m)find(x >= m,1,'first'),featVector));
                    % plot(x(rMatch),f(rMatch),'o','MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i})
                    % rPlot = (x>=min(featVector) & x<=max(featVector));
                    % if sum(rPlot) > 0
                    %     linePlots{i} = plot(x(rPlot),f(rPlot),'color',colors{i},'LineWidth',2);
                    % else
                    %     linePlots{i} = plot(ones(2,1)*x(rMatch(1)),[0,f(rMatch(1))],'color',colors{i},'LineWidth',2);
                    % end
                end
                % Add a legend if necessary
                if opi==1
                    legend([linePlots{:}],groupNames,'interpreter','none')
                end
                % % Add dots:
                % for i = 1:numGroups
                %     featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                %     [f, x] = ksdensity(featVector);
                %
                % end
                ylabel('Probability density')
            case 'hist'
                for i = 1:numGroups
                    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                    histogram(featVector,'BinMethod','auto','FaceColor',colors{i},'Normalization','probability')
                end
                if opi==1
                    legend(groupNames,'interpreter','none')
                end
                ylabel('Probability')
            end
            xlabel(sprintf('[%u] %s (%s=%4.2f%s)',Operations(op_ind).ID,Operations(op_ind).Name,...
                                cfnName,testStat(op_ind),cfnUnit),'interpreter','none')
        end
    end
end

%-------------------------------------------------------------------------------
% Dependence of top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'cluster'))

    % 1. Get pairwise similarity matrix
    op_ind = ifeat(1:numTopFeatures); % plot these operations indices

    % (if it exists already, use that; otherwise compute it on the fly)

    fileVarsStruct = whos('-file',whatDataFile);
    fileVars = {fileVarsStruct.name};
    if ismember('op_clust',fileVars)
        tmp = load(whatDataFile,'op_clust');
        clustStruct = tmp.op_clust; clear tmp
    else
        % Doesn't exist in the hctsa dataset:
        clustStruct = struct();
    end
    if isfield(clustStruct,'Dij') && ~isempty(clustStruct.Dij)
        % pairwise distances already computed, stored in the HCTSA .mat file
        fprintf(1,'Loaded %s distances from %s\n',clustStruct.distanceMetric,whatDataFile);
        Dij = squareform(clustStruct.Dij);
        Dij = Dij(op_ind,op_ind);
        distanceMetric = clustStruct.distanceMetric;
    else
        % Compute correlations on the fly
        Dij = BF_pdist(TS_DataMat(:,op_ind)','abscorr');
        distanceMetric = 'abscorr';
    end
    makeLabel = @(x) sprintf('[%u] %s (%4.2f%s)',Operations(x).ID,Operations(x).Name,...
                        testStat(x),cfnUnit);
    objectLabels = arrayfun(@(x)makeLabel(x),op_ind,'UniformOutput',0);
    BF_ClusterDown(Dij,floor(numTopFeatures/5),'whatDistance',distanceMetric,...
                        'objectLabels',objectLabels);
    title(sprintf('Dependencies between %u top features',numTopFeatures))
end

% Don't display lots of crap to screen unless the user wants it:
if nargout == 0
    clear ifeat testStat testStat_rand
end

%-------------------------------------------------------------------------------
% Functions
%-------------------------------------------------------------------------------
function tStat = fn_tStat(d1,d2)
    % Return test statistic from a 2-sample Welch's t-test
    [~,~,~,stats] = ttest2(d1,d2,'Vartype','unequal');
    tStat = stats.tstat;
end
%-------------------------------------------------------------------------------
function testStat = giveMeStats(dataMatrix,groupLabels)
    % Return test statistic for each operation
    testStat = zeros(numOps,1);
    for k = 1:numOps
        try
            testStat(k) = fn_testStat(dataMatrix(:,k),groupLabels,dataMatrix(:,k),groupLabels);
        catch
            testStat(k) = NaN;
        end
    end
end
%-------------------------------------------------------------------------------

end
