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
    whatTestStat = 'linclass';
    fprintf(1,'Using ''%s'' by default\n', whatTestStat);
end
if nargin < 3
    doNull = 1; % compute an empirical null distribution by randomizing class labels
end

% Use an inputParser to control plotting options:
inputP = inputParser;
default_whatPlots = {'histogram','distributions','cluster'};
check_whatPlots = @(x) iscell(x) || ischar(x);
addParameter(inputP,'whatPlots',default_whatPlots,check_whatPlots);
default_numTopFeatures = 25;
addParameter(inputP,'numTopFeatures',default_numTopFeatures,@isnumeric);
parse(inputP,varargin{:});

whatPlots = inputP.Results.whatPlots;
numTopFeatures = inputP.Results.numTopFeatures;
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
switch whatTestStat
case {'linear','linclass'}
    fprintf(1,'Using a linear classifier\n');
    cfnName = 'linear classifier';
    cfnUnit = '%';
    fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                    mean(yTest == classify(Xtest,XTrain,yTrain,'linear'))*100;
    chanceLine = 100/length(unique(timeSeriesGroup));
case 'diaglinear'
    fprintf(1,'A Naive Bayes classifier\n');
    cfnName = 'naive bayes classifier';
    cfnUnit = '%';
    fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                    mean(yTest == classify(Xtest,XTrain,yTrain,'diaglinear'))*100;
    chanceLine = 100/length(unique(timeSeriesGroup));
case {'svm','svmlinear'}
    fprintf(1,'A linear support vector machine\n');
    cfnName = 'SVM classifier';
    cfnUnit = '%';
    fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                    mean(yTest == predict(fitcsvm(XTrain,yTrain, ...
                            'KernelFunction','linear'),Xtest))*100;
    chanceLine = 100/length(unique(timeSeriesGroup));

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
[testStat_sort, ifeat] = sort(testStat,'descend');

isNaN = isnan(testStat_sort);
testStat_sort = testStat_sort(~isNaN);
ifeat = ifeat(~isNaN);
topN = min(20,length(Operations));
for i = 1:topN
    fprintf(1,'[%u] %s (%s) -- %4.2f%%\n',Operations(ifeat(i)).ID, ...
            Operations(ifeat(i)).Name,Operations(ifeat(i)).Keywords,testStat_sort(i));
end

% --------------------------------------------------------------------------
%% Plot outputs
% --------------------------------------------------------------------------

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
        fprintf(1,' %u %s statistics computed in %s.\n',numOps*numRepeats,...
                                        cfnName,BF_thetime(toc(timer)));
    end

    f = figure('color','w'); hold on
    colors = BF_getcmap('spectral',5,1);
    if ~doNull
        h_real = histogram(testStat,'Normalization','probability',...
                    'BinMethod','auto','FaceColor',colors{5},'FaceAlpha',0);
        maxH = max(h_real.Values);
    else
        % Plot both real distribution, and null distribution:
        numBins = 20;
        allTestStat = [testStat(:);testStat_rand(:)];
        minMax = [min(allTestStat),max(allTestStat)];
        histEdges = linspace(minMax(1),minMax(2),numBins+1);
        h_null = histogram(testStat_rand(:),histEdges,'Normalization','probability','FaceColor',colors{1});
        h_real = histogram(testStat,histEdges,'Normalization','probability','FaceColor',colors{5});
        maxH = max([max(h_real.Values),max(h_null.Values)]);
        legend('null','real')
        plot(mean(testStat_rand(:))*ones(2,1),[0,maxH],'--','color',colors{1},'LineWidth',2)
    end

    % Add chance line
    plot(chanceLine*ones(2,1),[0,maxH],'--k')
    % Add mean of real distribution
    plot(mean(testStat)*ones(2,1),[0,maxH],'--','color',colors{5},'LineWidth',2)
    % Labels:
    xlabel(sprintf('Individual %s across %u features',cfnName,numOps))
    ylabel('Probability')
end

%-------------------------------------------------------------------------------
% Distributions across classes for top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'distributions'))
    subPerFig = 5; % subplots per figure
    ks_or_hist = 'ks'; % view as either histograms or kernel-smoothed distributions
    colors = BF_getcmap('set1',max(numGroups,5),1);
    if numGroups==2, colors = colors([2,4]); end

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
                for i = 1:numGroups
                    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                    [f, x] = ksdensity(featVector);
                    plot(x,f,'color',colors{i},'LineWidth',2);
                end
                % Add a legend if necessary
                if opi==1
                    legend(groupNames)
                end
                % Add dots:
                for i = 1:numGroups
                    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                    [f, x] = ksdensity(featVector);
                    r = (arrayfun(@(m)find(x >= m,1,'first'),featVector));
                    plot(x(r),f(r),'o','MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i})
                end
                ylabel('Probability density')
            case 'hist'
                for i = 1:numGroups
                    featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
                    histogram(featVector,'BinMethod','auto','FaceColor',colors{i},'Normalization','probability')
                end
                if opi==1
                    legend(groupNames)
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
    else
        % Compute correlations on the fly
        Dij = BF_pdist(TS_DataMat(:,op_ind)','abscorr');
    end
    makeLabel = @(x) sprintf('[%u] %s (%4.2f)',Operations(x).ID,Operations(x).Name,...
                        testStat(x));
    objectLabels = arrayfun(@(x)makeLabel(x),op_ind,'UniformOutput',0);
    BF_ClusterDown(Dij,floor(numTopFeatures/5),'whatDistance','general',...
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
