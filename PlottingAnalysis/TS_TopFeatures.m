function [ifeat,testStat,testStat_rand,featureClassifier] = TS_TopFeatures(whatData,whatTestStat,varargin)
% TS_TopFeatures  Top individual features for discriminating labeled time series
%
% This function compares each feature in an hctsa dataset individually for its
% ability to separate the labeled classes of time series according to a given
% test statistic.
%
% Can also compare this performance to a set of randomized null features to
% evaluate the statistical significance of the result (pooled permutation testing).
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData, default: 'raw')
% whatTestStat, the test statistic to quantify the goodness of each feature
%               (e.g., 'fast_linear', 'tstat', 'svm', 'linear', 'diaglinear',
%                or others supported by GiveMeCfn)
%
%---OPTIONAL extra inputs:
%
% 'whatPlots', can specify what output plots to produce (cell of strings), e.g.,
%               specify {'histogram','distributions','cluster','datamatrix'} to
%               produce all four possible output plots (this is the default).
% 'numTopFeatures', can specify the number of top features to analyze, both in
%                   terms of the list of outputs, the histogram plots, and the
%                   cluster plot.
% 'numFeaturesDistr', can set a custom number of distributions to display (can
%                   set this lower to avoid producing large numbers of figures).
% 'numFolds', the number of folds used in cross-validation
% 'numNulls', the number of shuffled nulls to generate (e.g., 10 shuffles pools
%               shuffles for all M features, for a total of 10*M elements in the
%               null distribution) [default: 0]
% 'classifierFilename', .mat file to save the classifier to (not saved if empty).
%
%---EXAMPLE USAGE:
%
% TS_TopFeatures('norm','tstat','whatPlots',{'histogram','distributions',...
%           'cluster','datamatrix'},'numTopFeatures',40,'numFeaturesDistr',10);
%
%---OUTPUTS:
% ifeat, the ordering of operations by their performance
% testStat, the test statistic (whatTestStat) for each operation
% testStat_rand, test statistics making up the null distributions

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'raw';
end
if nargin < 2 || isempty(whatTestStat)
    whatTestStat = 'fast_linear'; % Way faster than proper prediction models
    fprintf(1,'Using ''%s'' test statistic by default\n', whatTestStat);
end

% Use an inputParser to control additional options as parameters:
inputP = inputParser;
% whatPlots
default_whatPlots = {'histogram','distributions','cluster'}; % 'datamatrix'
check_whatPlots = @(x) iscell(x) || ischar(x);
addParameter(inputP,'whatPlots',default_whatPlots,check_whatPlots);
% numTopFeatures
default_numTopFeatures = 40;
addParameter(inputP,'numTopFeatures',default_numTopFeatures,@isnumeric);
% numFeaturesDistr
default_numFeaturesDistr = 16;
addParameter(inputP,'numFeaturesDistr',default_numFeaturesDistr,@isnumeric);
% numNulls
default_numNulls = 0; % by default, don't compute an empirical null distribution
                      % by randomizing class labels
addParameter(inputP,'numNulls',default_numNulls,@isnumeric);
% numFolds
default_numFolds = [];
check_numFolds = @(x) isnumeric(x);
addParameter(inputP,'numFolds',default_numFolds,check_numFolds);
% classifierFilename
default_classifierFilename = '';
check_classifierFilename = @(x) ischar(x);
addParameter(inputP,'classifierFilename',default_classifierFilename,check_classifierFilename);

parse(inputP,varargin{:});

whatPlots = inputP.Results.whatPlots;
if ischar(whatPlots)
    whatPlots = {whatPlots};
end
numTopFeatures = inputP.Results.numTopFeatures;
numFeaturesDistr = inputP.Results.numFeaturesDistr;
numNulls = inputP.Results.numNulls;
numFolds = inputP.Results.numFolds;
classifierFilename = inputP.Results.classifierFilename;
clear('inputP');

%-------------------------------------------------------------------------------
%% Load the data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
numOps = height(Operations);
numTopFeatures = min(numTopFeatures,numOps);

%-------------------------------------------------------------------------------
%% Check that grouping information exists:
%-------------------------------------------------------------------------------
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
numClasses = max(TimeSeries.Group); % Assuming classes labeled with integers starting at 1
groupNames = TS_GetFromData(whatData,'groupNames');
if isempty(groupNames)
    error('No group label info in the data source');
end

%-------------------------------------------------------------------------------
% Fit the model
%-------------------------------------------------------------------------------
if isempty(numFolds) || numFolds==0
    % Use a heuristic to set a default number of folds
    numFolds = HowManyFolds(TimeSeries.Group,numClasses);
end

%-------------------------------------------------------------------------------
%% Define the train/test classification rate function, fn_testStat
%-------------------------------------------------------------------------------
% Also the chanceLine -- where you'd expect by chance (for equiprobable groups...)

switch whatTestStat
    case {'linear','linclass','fast_linear','diaglinear','svm','svm_linear'}
        % Set up the loss function for a classifier-based metric

        % (first check for possible class imbalance):
        classNumbers = arrayfun(@(x)sum(TimeSeries.Group==x),1:numClasses);
        isBalanced = all(classNumbers==classNumbers(1));
        if isBalanced
            fn_testStat = GiveMeFunctionHandle(whatTestStat,numClasses,'acc',false);
            fprintf(1,'Using total classification accuracy as output measure\n');
        else
            fn_testStat = GiveMeFunctionHandle(whatTestStat,numClasses,'balancedAcc',true);
            fprintf(1,'Due to class imbalance, using balanced classification accuracy as output measure\n');
        end
        chanceLine = 100/numClasses;
    case {'ustat','ranksum'}
        fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                                fn_uStat(XTrain(yTrain==1),XTrain(yTrain==2),false);
        chanceLine = NaN;
    case {'ustatExact','ranksumExact'}
        fn_testStat = @(XTrain,yTrain,Xtest,yTest) ...
                                fn_uStat(XTrain(yTrain==1),XTrain(yTrain==2),true);
        chanceLine = NaN;
    case {'ttest','tstat'}
        fn_testStat = @(XTrain,yTrain,Xtest,yTest,numFolds) ...
                                fn_tStat(XTrain(yTrain==1),XTrain(yTrain==2));
        chanceLine = 0; % chance-level t statistic is zero
    otherwise
        error('Unknown test statistics, ''%s''',whatTestStat);
end

% Now get information about the statistic and its units to display:
switch whatTestStat
    case {'linear','linclass','fast_linear'}
        testStatText = 'linear classifier';
        statUnit = '%';
    case 'diaglinear'
        testStatText = 'Naive bayes classifier';
        statUnit = '%';
    case {'svm','svm_linear'}
        testStatText = 'linear SVM classifier';
        statUnit = '%';
    case {'ustat','ranksum'}
        testStatText = 'Mann-Whitney approx p-value';
        statUnit = ' (log10(p))';
    case {'ustatExact','ranksumExact'}
        testStatText = 'Mann-Whitney exact p-value';
        statUnit = ' (log10(p))';
    case {'ttest','tstat'}
        testStatText = 'Welch''s t-stat';
        statUnit = ' (log10(p))';
        if numClasses > 2
            error('Cannot use t-test as test statistic with more than two groups :/');
        end
    otherwise
        error('Unknown method ''%s''',whatTestStat)
end

%-------------------------------------------------------------------------------
%% Loop over all features
%-------------------------------------------------------------------------------
% Use the same data for training and testing:
fprintf(1,'Comparing the (in-sample) performance of %u operations for %u classes using a %s...\n',...
                                height(Operations),numClasses,testStatText);
timer = tic;
testStat = giveMeStats(TS_DataMat,TimeSeries.Group,true);
fprintf(1,' Done in %s.\n',BF_TheTime(toc(timer)));

if all(isnan(testStat))
    error('Error computing test statistics for %s (maybe there are missing data?)',...
                testStatText);
end

%-------------------------------------------------------------------------------
% Give mean and that expected from random classifier (there may be a little overfitting)
if ~isnan(chanceLine)
    fprintf(1,['Mean %s performance across %u features = %4.2f%s\n' ...
            '(Random guessing for %u equiprobable classes = %4.2f%s)\n'], ...
        testStatText,numOps,nanmean(testStat),statUnit,numClasses,chanceLine,statUnit);
else
    fprintf(1,'Mean %s performance across %u features = %4.2f%s\n',...
        testStatText,numOps,nanmean(testStat),statUnit);
end

%-------------------------------------------------------------------------------
%% Display information about the top features (numTopFeatures)
%-------------------------------------------------------------------------------
[testStat_sort,ifeat] = sort(testStat,'descend'); % bigger is better

isNaN = isnan(testStat_sort);
testStat_sort = testStat_sort(~isNaN);
ifeat = ifeat(~isNaN);

% List the top features:
for i = 1:numTopFeatures
    ind = ifeat(i);
    fprintf(1,'[%u] %s (%s) -- %4.2f%s\n',Operations.ID(ind),...
            Operations.Name{ind},Operations.Keywords{ind},...
            testStat_sort(i),statUnit);
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%% Plot outputs
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Histogram of distribution of test statistics for labeled and null data
%-------------------------------------------------------------------------------
if ismember('histogram',whatPlots)
    % Plot the distribution of test statistics across all features:

    %-------------------------------------------------------------------------------
    %% Compute null distribution
    %-------------------------------------------------------------------------------
    testStat_rand = zeros(numOps,numNulls);
    if numNulls > 0
        fprintf(1,'Now for %u nulls... ',numNulls);
        tic
        for j = 1:numNulls
            if j < numNulls
                fprintf(1,'%u,',j);
            else
                fprintf(1,'%u.',j);
            end
            % Shuffle labels:
            groupLabels = TimeSeries.Group(randperm(height(TimeSeries)));
            testStat_rand(:,j) = giveMeStats(TS_DataMat,groupLabels,false);
        end
        fprintf(1,'\n%u %s statistics computed in %s.\n',numOps*numNulls,...
                                        testStatText,BF_TheTime(toc(timer)));

        % Pool nulls to estimate p-values
        pooledNulls = testStat_rand(:);
        pvals = arrayfun(@(x)mean(testStat(x) < pooledNulls),1:length(testStat));
        % FDR-corrected q-values:
        FDR_qvals = mafdr(pvals,'BHFDR','true');
        fprintf(1,'Estimating FDR-corrected p-values across all features by pooling across %u nulls\n',numNulls);
        fprintf(1,'(Given strong dependences across %u features, will produce conservative p-values)\n',numOps);
        % Give summary:
        sigThreshold = 0.05;
        if any(FDR_qvals < sigThreshold)
            fprintf(1,['%u/%u features show better performance using %s than the null distribution' ...
                        '\nat the magical 0.05 threshold (FDR corrected)\n'],...
                            sum(FDR_qvals < 0.05),length(FDR_qvals),testStatText);
        else
            fprintf(1,['Tough day at the office, hey? No features show statistically better performance than ' ...
                    'the null distribution at a FDR of 0.05.\nDon''t you go p-hacking now, will you?\n']);
        end
    end

    %---------------------------------------------------------------------------
    % Plot histogram
    f = figure('color','w'); hold on
    colors = BF_GetColorMap('spectral',5,1);
    if numNulls == 0
        % Just plot the real distribution of test statistics across all features
        h_real = histogram(testStat,'Normalization','probability',...
                    'BinMethod','auto','FaceColor',colors{5},'EdgeColor','k','FaceAlpha',0);
        maxH = max(h_real.Values);
    else
        % Plot both real distribution and null distribution:
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
    l_chance = plot(chanceLine*ones(2,1),[0,maxH],'--','color',colors{1});

    % Add mean of real distribution:
    l_mean = plot(nanmean(testStat)*ones(2,1),[0,maxH],'--','color',colors{5},'LineWidth',2);

    % Labels:
    xlabel(sprintf('Individual %s performance across %u features',testStatText,numOps))
    ylabel('Probability')

    % Legend:
    if numNulls > 0
        legend([h_null,h_real,l_chance,l_meannull,l_mean],'null','real','chance','null mean','mean')
        title(sprintf('%u features significantly informative of groups (FDR-corrected at 0.05)',sum(FDR_qvals < 0.05)))
    else
        legend([h_real,l_chance,l_mean],{'real','chance','mean'});
    end
else
    testStat_rand = [];
end

%-------------------------------------------------------------------------------
% Distributions across classes for top features
%-------------------------------------------------------------------------------
if ismember('distributions',whatPlots)
    subPerFig = 16; % subplots per figure

    % Set the colors to be assigned to groups:
    colors = GiveMeColors(numClasses);

    % Space the figures out properly:
    numFigs = ceil(numFeaturesDistr/subPerFig);

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
            r = r(r<=numFeaturesDistr);
        end
        featHere = ifeat(r); % features to plot on this figure
        % Make the figure
        f = figure('color','w');
        f.Position(3:4) = [1353, 857];
        % Loop through features
        for opi = 1:length(featHere)
            subplot(ceil(length(featHere)/4),4,opi);
            TS_SingleFeature(data,Operations.ID(featHere(opi)),true,false,...
                                                testStat(featHere(opi)),false);
        end
    end
end

%-------------------------------------------------------------------------------
% Data matrix containing top features
%-------------------------------------------------------------------------------
if ismember('datamatrix',whatPlots)
    featInd = ifeat(1:numTopFeatures);
    ixFeat = BF_ClusterReorder(TS_DataMat(:,featInd)','corr','average');
    dataLocal = struct('TS_DataMat',BF_NormalizeMatrix(TS_DataMat(:,featInd(ixFeat)),'maxmin'),...
                    'TimeSeries',TimeSeries,...
                    'Operations',Operations(featInd(ixFeat),:));
    TS_PlotDataMatrix(dataLocal,'colorGroups',true,'groupReorder',true);
end

%-------------------------------------------------------------------------------
% Inter-dependence of top features
%-------------------------------------------------------------------------------
if ismember('cluster',whatPlots)

    % 1. Get pairwise similarity matrix
    op_ind = ifeat(1:numTopFeatures); % plot these operations indices

    % (if it exists already, use that; otherwise compute it on the fly)
    clustStruct = TS_GetFromData(whatData,'op_clust');
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
    makeLabel = @(x) sprintf('[%u] %s (%4.2f%s)',Operations.ID(x),Operations.Name{x},...
                        testStat(x),statUnit);
    objectLabels = arrayfun(@(x)makeLabel(x),op_ind,'UniformOutput',0);
    clusterThreshold = 0.2; % threshold at which split into clusters
    [~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',objectLabels);
    title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                            numTopFeatures,length(cluster_Groupi)))
end

%-------------------------------------------------------------------------------
% Save a classifier (for the single best feature) to file
%-------------------------------------------------------------------------------
featureClassifier = struct();
if ~isempty(classifierFilename)
    theTopFeatureIndex = ifeat(1);
    topFeatureValues = TS_DataMat(:,theTopFeatureIndex);
    numFoldsNow = 0;

    % Train a classification model on the top feature:
    [bestTestStat,bestMdl,whatTestStat] = fn_testStat(topFeatureValues,...
                TimeSeries.Group,topFeatureValues,TimeSeries.Group,numFoldsNow);

    % Prepare the featureClassifier structure:
    featureClassifier.Operation.ID = Operations.ID(theTopFeatureIndex); % Sorted list of top operations
    featureClassifier.Operation.Name = Operations.Name{theTopFeatureIndex}; % Sorted list of top operations
    featureClassifier.Mdl = bestMdl; % sorted feature models
    if numFolds > 0
        % From sorted accuracy of models
        featureClassifier.CVAccuracy = testStat_sort(1);
    else
        % CV was never done
        featureClassifier.CVAccuracy = NaN;
    end
    featureClassifier.Accuracy = bestTestStat;
    featureClassifier.whatTestStat = whatTestStat;
    featureClassifier.normalizationInfo = TS_GetFromData(whatData,'normalizationInfo');
    classes = groupNames;

    % Save to file
    if exist(classifierFilename,'file')~=0
        save(classifierFilename,'featureClassifier','classes','-append');
        fprintf(1,'Appended individual-feature classifiers to %s\n',classifierFilename);
    else
        save(classifierFilename,'featureClassifier','classes','-v7.3');
        fprintf(1,'Saved individual-feature classifiers to %s\n',classifierFilename);
    end
end

%-------------------------------------------------------------------------------
% Don't display crap to screen unless the user wants it:
if nargout == 0
    clear('ifeat','testStat','testStat_rand','featureClassifier');
end

%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
function [uStatP,Mdl] = fn_uStat(d1,d2,doExact)
    % Return test statistic from Mann-Whitney U-test
    if doExact
        [p,~,stats] = ranksum(d1,d2,'method','exact');
    else
        [p,~,stats] = ranksum(d1,d2,'method','approx');
    end
    uStatP = -log10(p);
    % uStat = stats.ranksum;
    Mdl = stats;
end
%-------------------------------------------------------------------------------
function [tStat,Mdl] = fn_tStat(d1,d2)
    % Return test statistic from a 2-sample Welch's t-test
    [~,~,~,stats] = ttest2(d1,d2,'Vartype','unequal');
    tStat = stats.tstat;
    Mdl = stats;
end
%-------------------------------------------------------------------------------
function [testStat,Mdl] = giveMeStats(dataMatrix,groupLabels,beVerbose)
    % Return test statistic for each operation
    testStat = zeros(numOps,1);
    Mdl = cell(numOps,1);
    loopTimer = tic;
    for k = 1:numOps
        try
            if nargout == 2
                % This is slower for the fast_linear classifier (but returns a model)
                [testStat(k),Mdl{k}] = fn_testStat(dataMatrix(:,k),groupLabels,dataMatrix(:,k),groupLabels);
            else
                testStat(k) = fn_testStat(dataMatrix(:,k),groupLabels,dataMatrix(:,k),groupLabels);
            end
        catch
            % keyboard
            fprintf('Could not return model for operation %u',k);
        end
        % Give estimate of time remaining:
        if beVerbose && k==100
            fprintf(1,'(should take approx %s to compute for all %u features)\n',...
                            BF_TheTime(toc(loopTimer)/100*(numOps)),numOps);
        end
    end
end
%-------------------------------------------------------------------------------

end
