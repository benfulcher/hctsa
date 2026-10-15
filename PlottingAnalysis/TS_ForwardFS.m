function [ifeat,testStat,trainErr,testErr,testClass] = TS_ForwardFS(whatData,iTrain,criterion,cvFolds,numFeatSelect,annotateParams)
% TS_ForwardFS  Greedy forward feature selection
%
% Uses the sequentialfs function from Matlab's Statistics Toolbox to select
% features one at a time, at each step adding the feature that maximizes the
% cross-validated classification rate of the training data (using the specified
% training indices).
% After selecting the features, the learned classification rule is applied to the
% training and test sets (all data not in the training set) to get training and
% test classification errors for the feature set after each feature is added.
% The first two selected features are plotted against each other.
%
% Typical usage uses 'fast-linear' for the criterion (linear classification rates),
% or 'diaglinear' when there are highly correlated features in the data matrix.
%
% NOTE: A training portion should be specified. If there
%       is no obvious training portion of your data, this should be specified
%       at random (e.g., stratified to represent groups in your data) and
%       repeated to get an estimate of the out-of-sample classification that is
%       not dependent on the particular data partition.
%       If iTrain is empty, all (labeled) data is used for training, and there is
%       no test set.
%
%---INPUTS:
% whatData: the data to load in (cf. TS_LoadData)
% iTrain: indices of data to train on (rows of the data matrix; default: all data)
% criterion: what criterion on which to evaluate the quality of a feature set
%           (a classifier, cf. GiveMeCfn)
% cvFolds: the number of cross-validation folds to use in the feature selection
%           process
% numFeatSelect: the total number of features to select
% annotateParams: annotation parameters for the 2-d space
%
%---OUTPUTS:
% ifeat: indices of features selected, in the order selected.
% testStat: cross-validated classification rate (%) of the training data, for the
%           feature set after each feature is added.
% trainErr: training errors (%) for selected features, after each feature is added.
% testErr: test errors (%) for selected features, after each feature is added
%           (empty if there is no test set).
% testClass: classification of the test data using all selected features (empty if
%           there is no test set).

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
end

if nargin < 2
    iTrain = [];
end

if nargin < 3 || isempty(criterion)
    criterion = 'fast-linear';
    fprintf(1,'Default test statistic: Using fast, in-sample linear classification rate\n');
end

if nargin < 4 || isempty(cvFolds)
    cvFolds = 5;
    fprintf(1,'Default: 5-fold cross-validation\n');
end

if nargin < 5 || isempty(numFeatSelect)
    numFeatSelect = 2; % Stop after two features are selected
end

if nargin < 6
    annotateParams = struct('n',6,'textAnnotation','Name');
end

%-------------------------------------------------------------------------------
%% Load the data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

% Training indices refer to rows of the loaded data; keep track of them by ID
% through removing unlabeled time series:
if isempty(iTrain)
    iTrain = 1:height(TimeSeries);
end
trainIDs = TimeSeries.ID(iTrain);
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
isTrain = ismember(TimeSeries.ID,trainIDs);
isTest = ~isTrain;
hasTest = any(isTest);

XTrain = TS_DataMat(isTrain,:);
yTrain = TimeSeries.Group(isTrain);

% ------------------------------------------------------------------------------
%% Set up the classification function
% ------------------------------------------------------------------------------
% Parameters for evaluating classification rates (in-sample or cross-validated):
cfnParams = GiveMeDefaultClassificationParams(TimeSeries(isTrain,:),[],false);
cfnParams.whatClassifier = criterion;
cfnParams = UpdateClassifierText(cfnParams);
cfnParams.numFolds = 0;
cfnParams.computePerFold = false;

% sequentialfs needs a loss, summed over observations, to minimize:
selectParams = cfnParams;
if cfnParams.doReweight
    selectParams.whatLoss = 'balancedLoss';
    fprintf(1,'Unbalanced classes: using a balanced accuracy measure to select features (& using reweighting)...\n');
else
    selectParams.whatLoss = 'sumLoss';
end
classify_fn = GiveMeFunctionHandle(selectParams);

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
%% Do the feature selection
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
FS_timer = tic;
fprintf(1,'Performing greedy forward feature selection using ''%s''...\n',criterion);

opts = statset('display','iter');
[fs,history] = sequentialfs(classify_fn,XTrain,yTrain,...
                    'cv',cvFolds,'options',opts,'nfeatures',numFeatSelect);
fprintf(1,'Selected %u features\n',sum(fs));

% Finished selecting features!
fprintf(1,'Feature selection to %u features completed in %s.\n',...
            numFeatSelect,BF_TheTime(toc(FS_timer)))
clear('FS_timer');

% Features in the order that they were selected:
numSelected = size(history.In,1);
ifeat = zeros(1,numSelected);
ifeat(1) = find(history.In(1,:));
for i = 2:numSelected
    ifeat(i) = find(history.In(i,:) & ~history.In(i-1,:));
end

%-------------------------------------------------------------------------------
%% Evaluate the selected feature sets
%-------------------------------------------------------------------------------
% Cross-validated classification rate on the training data:
cvParams = cfnParams;
cvParams.numFolds = cvFolds;

testStat = zeros(numSelected,1);
trainErr = zeros(numSelected,1);
testErr = zeros(numSelected,1);
for i = 1:numSelected
    XTrain_i = XTrain(:,ifeat(1:i));
    testStat(i) = GiveMeCfn(XTrain_i,yTrain,[],[],cvParams);
    trainErr(i) = 100 - GiveMeCfn(XTrain_i,yTrain,XTrain_i,yTrain,cfnParams);
    if hasTest
        testErr(i) = 100 - GiveMeCfn(XTrain_i,yTrain,TS_DataMat(isTest,ifeat(1:i)),...
                                    TimeSeries.Group(isTest),cfnParams);
    end
end

% Classify the test data using the full set of selected features:
if hasTest
    [~,Mdl] = GiveMeCfn(XTrain(:,ifeat),yTrain,[],[],cfnParams);
    testClass = predict(Mdl,TS_DataMat(isTest,ifeat));
else
    testErr = [];
    testClass = [];
end

%-------------------------------------------------------------------------------
% Plot in a 2-D space
%-------------------------------------------------------------------------------
if numSelected > 1
    % Set feature labels:
    featureLabels = cell(2,1);
    for i = 1:2
        featureLabels{i} = sprintf('%s (%.2f%%)',Operations.Name{ifeat(i)},testStat(i));
    end
    TS_Plot2d(TS_DataMat(:,ifeat(1:2)),TimeSeries,featureLabels,annotateParams);
end

end
