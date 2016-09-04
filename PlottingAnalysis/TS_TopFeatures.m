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
%               (e.g., 'fast_linear', 'tstat', 'svm', 'linear', 'diaglinear',
%                or others supported by GiveMeFunctionHandle)
% doNull, (boolean) whether to compute a null distribution using permutations
%           of the class labels
%
% ***Additional plotting options***:
% 'whatPlots', can specify what output plots to produce (cell of strings), e.g.,
%               specify {'histogram','distributions','cluster','datamatrix'} to
%               produce all four possible output plots (this is the default).
% 'numTopFeatures', can specify the number of top features to analyze, both in
%                   terms of the list of outputs, the histogram plots, and the
%                   cluster plot.
% 'numHistogramFeatures', can optionally also set a custom number of histograms
%                       to display (often want to set this lower to avoid producing
%                       large numbers of figures).
%
%---EXAMPLE USAGE:
%
% TS_TopFeatures('norm','tstat',1,'whatPlots',{'histogram','distributions',...
%           'cluster','datamatrix'},'numTopFeatures',40,'numHistogramFeatures',10);
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
    doNull = 0; % compute an empirical null distribution by randomizing class labels
end

% Use an inputParser to control additional plotting options as parameters:
inputP = inputParser;
default_whatPlots = {'histogram','distributions','cluster'}; % datamatrix
check_whatPlots = @(x) iscell(x) || ischar(x);
addParameter(inputP,'whatPlots',default_whatPlots,check_whatPlots);
default_numTopFeatures = 40;
addParameter(inputP,'numTopFeatures',default_numTopFeatures,@isnumeric);
default_numHistogramFeatures = 16;
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
numClasses = max(timeSeriesGroup); % Assuming classes labeled with integers starting at 1
groupNames = TS_GetFromData(whatDataFile,'groupNames');
if isempty(groupNames)
    error('No group label info in the data source');
end

% --------------------------------------------------------------------------
%% Define the train/test classification rate function, fn_testStat
% --------------------------------------------------------------------------
% Also the chanceLine -- where you'd expect by chance (for equiprobable groups...)
if ismember(whatTestStat,{'linear','linclass','fast_linear','diaglinear','svm','svm_linear'})

    %-------------------------------------------------------------------------------
    % Set up the loss function

    % (first check for possible class imbalance):
    classNumbers = arrayfun(@(x)sum(timeSeriesGroup==x),1:numClasses);
    isBalanced = all(classNumbers==classNumbers(1));
    if isBalanced
        fn_testStat = GiveMeFunctionHandle(whatTestStat,numClasses,'acc',0);
        fprintf(1,'Using overall classification accuracy as output measure\n');
    else
        fn_testStat = GiveMeFunctionHandle(whatTestStat,numClasses,'balancedAcc',1);
        fprintf(1,'Due to class imbalance, using balanced classification accuracy as output measure\n');
    end

    % fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                    % mean(yTest == classify(Xtest,XTrain,yTrain,'linear'))*100;
    chanceLine = 100/numClasses;
    cfnUnit = '%';
end

switch whatTestStat
case {'linear','linclass','fast_linear'}
    cfnName = 'linear classifier';
case 'diaglinear'
    cfnName = 'Naive bayes classifier';
case {'svm','svm_linear'}
    cfnName = 'linear SVM classifier';
case {'ttest','tstat'}
    cfnName = 'Welch''s t-stat';
    cfnUnit = '';
    if numClasses > 2
        error('Cannot use t-test as test statistic with more than two groups :/');
    end
    fn_testStat = @(XTrain,yTrain,Xtest,yTest) fn_tStat(XTrain(yTrain==1),XTrain(yTrain==2));
    chanceLine = 0; % by chance, t-stats averge to zero
otherwise
    error('Unknown method ''%s''',whatTestStat)
end

% --------------------------------------------------------------------------
%% Loop over all features
% --------------------------------------------------------------------------
% Use the same data for training and testing:
fprintf(1,'Comparing the (in-sample) performance of %u operations for %u classes using a %s...\n',...
                                length(Operations),numClasses,cfnName);
timer = tic;
testStat = giveMeStats(TS_DataMat,timeSeriesGroup,1);
fprintf(1,' Done in %s.\n',BF_thetime(toc(timer)));

% Give mean and that expected from random classifier (may be a little overfitting)
fprintf(1,['Mean %s performance across %u operations = %4.2f%s\n' ...
            '(Random guessing for %u equiprobable classes = %4.2f%s)\n'], ...
        cfnName,numOps,nanmean(testStat),cfnUnit,numClasses,chanceLine,cfnUnit);

% --------------------------------------------------------------------------
%% Display information about the top topN operations
% --------------------------------------------------------------------------
[testStat_sort, ifeat] = sort(testStat,'descend'); % bigger is better

isNaN = isnan(testStat_sort);
testStat_sort = testStat_sort(~isNaN);
ifeat = ifeat(~isNaN);

% List the top features:
topN = min(numTopFeatures,length(Operations));
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
    numNulls = 10;
    testStat_rand = zeros(numOps,numNulls);
    if doNull
        fprintf(1,'Now for %u nulls... ',numNulls);
        tic
        for j = 1:numNulls
            if j < numNulls
                fprintf(1,'%u,',j);
            else
                fprintf(1,'%u',j);
            end
            % Shuffle labels:
            groupLabels = timeSeriesGroup(randperm(length(timeSeriesGroup)));
            testStat_rand(:,j) = giveMeStats(TS_DataMat,groupLabels,0);
        end
        fprintf(1,'\n%u %s statistics computed in %s.\n',numOps*numNulls,...
                                        cfnName,BF_thetime(toc(timer)));

        % Pool nulls to estimate p-values
        pooledNulls = testStat_rand(:);
        pvals = arrayfun(@(x)mean(testStat(x) < pooledNulls),1:length(testStat));
        % FDR-corrected q-values:
        FDR_qvals = mafdr(pvals,'BHFDR','true');
        fprintf(1,'Estimating FDR-corrected p-values across all features by pooling across %u nulls\n',numNulls);
        fprintf(1,'(Given strong dependences across %u features, will produce conservative p-values)\n',numOps);
        % Give summary:
        if any(FDR_qvals < 0.05)
            fprintf(1,['%u/%u features show better performance using %s than the null distribution' ...
                        '\nat the magical 0.05 threshold (FDR corrected)\n'],...
                            sum(FDR_qvals<0.05),length(FDR_qvals),cfnName);
        else
            fprintf(1,['Tough day at the office, hey? No features show statistically better performance than ' ...
                    'the null distribution at a FDR of 0.05.\nDon''t you go p-hacking now, will you?\n']);
        end
    end

    f = figure('color','w'); hold on
    colors = BF_getcmap('spectral',5,1);
    if ~doNull
        h_real = histogram(testStat,'Normalization','probability',...
                    'BinMethod','auto','FaceColor',colors{5},'EdgeColor','k','FaceAlpha',0);
        maxH = max(h_real.Values);
    else
        % Plot both real distribution, and null distribution:
        numBins = 20;
        allTestStat = [testStat(:);testStat_rand(:)];
        minMax = [min(allTestStat),max(allTestStat)];
        histEdges = linspace(minMax(1),minMax(2),numBins+1);
        h_null = histogram(testStat_rand(:),histEdges,'Normalization','probability','FaceColor',colors{1});
        h_real = histogram(testStat,histEdges,'Normalization','probability',...
                                'FaceColor',colors{5},'EdgeColor','k');
        maxH = max([max(h_real.Values),max(h_null.Values)]);
        l_meannull = plot(nanmean(testStat_rand(:))*ones(2,1),[0,maxH],'--','color',colors{1},'LineWidth',2);
    end

    % Add chance line:
    l_chance = plot(chanceLine*ones(2,1),[0,maxH],'--k');

    % Add mean of real distribution:
    l_mean = plot(nanmean(testStat)*ones(2,1),[0,maxH],'--','color',colors{5},'LineWidth',2);

    % Labels:
    xlabel(sprintf('Individual %s across %u features',cfnName,numOps))
    ylabel('Probability')

    % Legend:
    if doNull
        legend([h_null,h_real,l_chance,l_meannull,l_mean],'null','real','chance','null mean','real mean')
    else
        legend([h_real,l_chance,l_mean],{'real','chance','real mean'});
    end
end

%-------------------------------------------------------------------------------
% Distributions across classes for top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'distributions'))
    subPerFig = 16; % subplots per figure
    % ks_or_hist = 'ks'; % 'ks'/'hist': view as either histograms or kernel-smoothed distributions

    % Set the colors to be assigned to groups:
    colors = GiveMeColors(numClasses);

    % Space the figures out properly:
    numFigs = ceil(numHistogramFeatures/subPerFig);

    % Make data structure for TS_SingleFeature
    data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);
    data.groupNames = groupNames;

    for figi = 1:numFigs
        if figi*subPerFig > length(ifeat)
            break % We've exceeded number of features
        end
        % Get the indices of features to plot
        r = ((figi-1)*subPerFig+1:figi*subPerFig);
        if figi==numFigs % filter down for last one
            r = r(r<=numHistogramFeatures);
        end
        featHere = ifeat(r); % features to plot on this figure
        % Make the figure
        f = figure('color','w');
        f.Position(3:4) = [1353, 857];
        % Loop through features
        for opi = 1:length(featHere)
            subplot(ceil(length(featHere)/4),4,opi);
            TS_SingleFeature(data,featHere(opi),1,0,testStat(featHere(opi)),0);
        end
    end

    %
    % for figi = 1:numFigs
    %     if figi*subPerFig > length(ifeat)
    %         break % We've exceeded number of features
    %     end
    %     r = ((figi-1)*subPerFig+1:figi*subPerFig);
    %     featHere = ifeat(r); % features to plot on this figure
    %
    %     f = figure('color','w');
    %     f.Position(3:4) = [588, 612]; % make longer
    %     for opi = 1:subPerFig
    %         op_ind = featHere(opi);
    %         ax = subplot(subPerFig,1,opi); hold on
    %
    %         switch ks_or_hist
    %         case 'ks'
    %             % Plot distributions first for the sake of the legend
    %             linePlots = cell(numClasses,1);
    %             for i = 1:numClasses
    %                 featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
    %                 [~,~,linePlots{i}] = BF_plot_ks(featVector,colors{i},0,2,20,1);
    %             end
    %             % Trim x-limits (with 2% overreach)
    %             ax.XLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
    %             ax.XLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));
    %             % Add a legend if necessary
    %             if opi==1
    %                 legend([linePlots{:}],groupNames,'interpreter','none')
    %             end
    %             ylabel('Probability density')
    %         case 'hist'
    %             for i = 1:numClasses
    %                 featVector = TS_DataMat((timeSeriesGroup==i),op_ind);
    %                 histogram(featVector,'BinMethod','auto','FaceColor',colors{i},'Normalization','probability')
    %             end
    %             if opi==1
    %                 legend(groupNames,'interpreter','none')
    %             end
    %             ylabel('Probability')
    %         end
    %
    %         % Annotate rectangles under the distributions
    %         BF_AnnotateRect(whatTestStat,TS_DataMat(:,op_ind),timeSeriesGroup,numClasses,colors,ax);
    %
    %         xlabel(sprintf('[%u] %s (%s=%4.2f%s)',Operations(op_ind).ID,Operations(op_ind).Name,...
    %                             cfnName,testStat(op_ind),cfnUnit),'interpreter','none')
    %     end
    % end
end

%-------------------------------------------------------------------------------
% Data matrix containing top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'datamatrix'))
    featInd = ifeat(1:numTopFeatures);
    ixFeat = BF_ClusterReorder(TS_DataMat(:,featInd)','corr','average');
    dataLocal = struct('TS_DataMat',BF_NormalizeMatrix(TS_DataMat(:,featInd(ixFeat)),'maxmin'),...
                    'TimeSeries',TimeSeries,...
                    'Operations',Operations(featInd(ixFeat)));
    TS_plot_DataMatrix(dataLocal,'colorGroups',1,'groupReorder',1);
end

%-------------------------------------------------------------------------------
% Inter-dependence of top features
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'cluster'))

    % 1. Get pairwise similarity matrix
    op_ind = ifeat(1:numTopFeatures); % plot these operations indices

    % (if it exists already, use that; otherwise compute it on the fly)
    clustStruct = TS_GetFromData(whatDataFile,'op_clust');
    if isempty(clustStruct) % doesn't exist
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
    clusterThreshold = 0.2; % threshold at which split into clusters
    [~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',objectLabels);
    title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                            numTopFeatures,length(cluster_Groupi)))
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
function testStat = giveMeStats(dataMatrix,groupLabels,beVerbose)
    % Return test statistic for each operation
    testStat = zeros(numOps,1);
    loopTimer = tic;
    for k = 1:numOps
        try
            testStat(k) = fn_testStat(dataMatrix(:,k),groupLabels,dataMatrix(:,k),groupLabels);
        catch
            testStat(k) = NaN;
        end
        % Give estimate of time remaining:
        if beVerbose && k==20
            fprintf(1,'(should take approx %s to compute for all %u features)\n',...
                            BF_thetime(toc(loopTimer)/20*(numOps)),numOps);
        end
    end
end
%-------------------------------------------------------------------------------

end
