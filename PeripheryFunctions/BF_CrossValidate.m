function [accuracyMean,accuracyStd,inSampleAccStd,allFoldData,cvOut] = BF_CrossValidate(X,y,cfnParams,spec,varargin)
% BF_CrossValidate   Repeated, stratified cross-validation with an optional
%                    transformation fitted within each training fold.
%
% For each repeat, the observations are split into cfnParams.numFolds stratified
% folds (cvpartition, using the current random seed). For each fold:
%   (i) if spec is given, a feature filter and normalization is fitted on the
%       training fold (BF_FitTransform) and applied to both folds
%       (BF_ApplyTransform), so the test fold plays no part in fitting it;
%   (ii) the classifier of cfnParams is trained on the training fold (GiveMeCfn),
%       with the class names, class prior and observation weights of all
%       observations, as in the built-in cross-validation of the fitc*
%       functions (these use only the labels);
%   (iii) it predicts the test fold.
% Performance (cfnParams.whatLoss) is summarized as in ComputeCVAccuracies, with
% outputs of the same form.
%
% Without spec (empty), the data are used as given: with the same partitions,
% the predictions (and so the performance) are those of the built-in
% cross-validation of the classifier ('KFold' or 'CVPartition' in fitcsvm,
% fitcecoc, ...) used by ComputeCVAccuracies, and from the same random seed the
% partitions are the same too. (Classification scores can differ slightly, within
% the convergence tolerance of iterative solvers such as fitcsvm's.)
%
%---USAGE:
% meanAcc = BF_CrossValidate(X,y,cfnParams); % no transformation
% meanAcc = BF_CrossValidate(X,y,cfnParams,struct('normFunction','mixedSigmoid'));
%
%---INPUTS:
% X, the data matrix (observations x features): raw (unnormalized) data if spec
%       is given.
% y, the class label of each observation (categorical).
% cfnParams, classification parameters (cf. GiveMeDefaultClassificationParams).
% spec, the normalization specification fitted within each training fold (cf.
%       BF_FitTransform): fields normFunction, filterOptions (only the column
%       proportion is used: apply the row filter beforehand), classVarFilter,
%       testNaN. Missing values of kept features in a test fold are replaced by
%       the median of the normalized training-fold values unless spec.testNaN
%       is 'keep' (default here: 'impute').
%       Default: [] (no transformation).
%
% (name-value options):
% 'Operations': a table with the ID and Name of each column of X (or a vector
%       of IDs), stored in the fitted transformations. Default: IDs 1, 2, ...
% 'partitions': the cvpartition of each repeat (a cell of numRepeats
%       cvpartition objects), e.g., cvOut.partitions of an earlier call.
%       Default: new stratified partitions, cvpartition(y,'KFold',numFolds).
% 'transforms': transformations fitted earlier on the training folds of these
%       partitions (cvOut.transforms of an earlier call with the same
%       'partitions'), used instead of fitting new ones. Only valid for a
%       transformation that does not depend on the labels (classVarFilter false).
% 'doSummarize': whether to average the outputs across repeats. Default: true.
%
%---OUTPUTS:
% accuracyMean, accuracyStd, inSampleAccStd, allFoldData: as for
%       ComputeCVAccuracies (one row per repeat, or the mean across repeats if
%       doSummarize; allFoldData is numFolds x 2 (training, test) x numRepeats).
% cvOut, a structure with fields:
%   partitions: the cvpartition of each repeat (empty if cfnParams.numFolds = 0)
%   transforms: the transformation fitted on each training fold
%               (numFolds x numRepeats cell; empty if no spec)
%   yPredict: the predicted label of each observation when it was in the test
%               fold (observations x numRepeats); in-sample predictions if
%               cfnParams.numFolds = 0
%   numFeatures: the number of features used in each fold (numFolds x numRepeats)
%   numImputed: the number of test-fold values imputed in each fold
%               (numFolds x numRepeats)

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
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 4
    spec = [];
end
inputP = inputParser;
addParameter(inputP,'Operations',[]);
addParameter(inputP,'partitions',{},@(x) iscell(x) || isa(x,'cvpartition'));
addParameter(inputP,'transforms',{},@iscell);
addParameter(inputP,'doSummarize',true,@(x) islogical(x) || isnumeric(x));
parse(inputP,varargin{:});
Operations = inputP.Results.Operations;
partitions = inputP.Results.partitions;
transforms = inputP.Results.transforms;
doSummarize = inputP.Results.doSummarize;
clear('inputP');
if isempty(Operations)
    Operations = (1:size(X,2))'; % column indices as feature IDs
end

y = y(:);
numObs = length(y);
if size(X,1)~=numObs
    error('X has %u rows, but there are %u labels',size(X,1),numObs);
end
numFolds = cfnParams.numFolds;
numRepeats = cfnParams.numRepeats;
if numFolds==0 && numRepeats~=1
    error('Without cross-validation (numFolds = 0), numRepeats must be 1');
end

% The transformation fitted within each training fold:
doTransform = ~isempty(spec);
if doTransform && (~isfield(spec,'testNaN') || isempty(spec.testNaN))
    spec.testNaN = 'impute';
end
reuseTransforms = ~isempty(transforms);
if reuseTransforms
    if ~doTransform
        error('''transforms'' were given without a normalization specification');
    end
    if isfield(spec,'classVarFilter') && ~isempty(spec.classVarFilter) && spec.classVarFilter
        error('Transformations that filter on class variance depend on the labels, so cannot be reused');
    end
    if isempty(partitions) && numFolds > 0
        error('''transforms'' can only be reused with the ''partitions'' they were fitted on');
    end
end

% Partitions:
if isa(partitions,'cvpartition')
    partitions = {partitions};
end
if numFolds > 0 && ~isempty(partitions)
    if length(partitions)~=numRepeats
        error('%u partitions given for %u repeats',length(partitions),numRepeats);
    end
    for r = 1:numRepeats
        if partitions{r}.NumObservations~=numObs || partitions{r}.NumTestSets~=numFolds
            error('Partition %u does not have %u folds of %u observations',r,numFolds,numObs);
        end
    end
end
if reuseTransforms && ~isequal(size(transforms),[max(numFolds,1),numRepeats])
    error('''transforms'' should be a %u x %u cell',max(numFolds,1),numRepeats);
end

% The classifier is trained on each training fold (GiveMeCfn without CV):
cfnTrain = cfnParams;
cfnTrain.numFolds = 0;

% Results are only computed per fold when there are folds:
perFold = cfnParams.computePerFold && (numFolds > 0);

%-------------------------------------------------------------------------------
%% Initialize
%-------------------------------------------------------------------------------
accuracyMean = zeros(numRepeats,1);
accuracyStd = nan(numRepeats,1);
inSampleAccStd = nan(numRepeats,2);
allFoldData = nan(numFolds,2,numRepeats); % folds, train/test, repeats
numSets = max(numFolds,1); % (a single training set without CV)
if isempty(partitions)
    partitions = cell(numRepeats,1);
end
if ~reuseTransforms
    transforms = cell(numSets,numRepeats);
    if ~doTransform
        transforms = {};
    end
end
yPredict = repmat(y,1,numRepeats);
numFeatures = zeros(numSets,numRepeats);
numImputed = zeros(numSets,numRepeats);

%-------------------------------------------------------------------------------
%% Cross-validate
%-------------------------------------------------------------------------------
for r = 1:numRepeats
    % Stratified partition (created here, so the random number stream is used in
    % the same order as by the built-in 'KFold' cross-validation of each repeat):
    if numFolds > 0 && isempty(partitions{r})
        partitions{r} = cvpartition(y,'KFold',numFolds);
    end
    accTrain = nan(numSets,1);
    accTest = nan(numSets,1);

    for k = 1:numSets
        % Training and test observations:
        if numFolds > 0
            isTrain = training(partitions{r},k);
            isTest = test(partitions{r},k);
        else
            isTrain = true(numObs,1);
            isTest = true(numObs,1);
        end

        % Transform the data, fitting only on the training fold:
        if doTransform
            if reuseTransforms
                T = transforms{k,r};
                XTrain = BF_ApplyTransform(X(isTrain,:),T,Operations,false);
            else
                [T,XTrain] = BF_FitTransform(X(isTrain,:),spec,Operations,y(isTrain));
                transforms{k,r} = T;
            end
            [XTest,isImputed] = BF_ApplyTransform(X(isTest,:),T,Operations);
            numImputed(k,r) = sum(isImputed(:));
        else
            XTrain = X(isTrain,:);
            XTest = X(isTest,:);
        end
        numFeatures(k,r) = size(XTrain,2);

        % Train on the training fold and predict the test fold:
        if numFolds > 0
            [~,Mdl] = GiveMeCfn(XTrain,y(isTrain),[],[],cfnTrain,false,y);
        else
            [~,Mdl] = GiveMeCfn(XTrain,y,[],[],cfnTrain);
        end
        yPredict(isTest,r) = predict(Mdl,XTest);

        if perFold
            accTest(k) = BF_LossFunction(y(isTest),yPredict(isTest,r),...
                                cfnParams.whatLoss,cfnParams.classLabels);
            accTrain(k) = BF_LossFunction(y(isTrain),predict(Mdl,XTrain),...
                                cfnParams.whatLoss,cfnParams.classLabels);
        end
    end

    if perFold
        allFoldData(:,1,r) = accTrain; % training folds
        allFoldData(:,2,r) = accTest; % test folds
        accuracyMean(r) = mean(accTest);
        accuracyStd(r) = std(accTest);
        inSampleAccStd(r,:) = [mean(accTrain),std(accTrain)];
        % Warning for training fold over-fitting:
        if ~cfnParams.suppressWarning && all(accTrain==100)
            warning('100% in-sample accuracy alert (just FYI: no need for any action!)')
        end
    else
        % Computed once from the predictions of all test folds (or in-sample):
        accuracyMean(r) = BF_LossFunction(y,yPredict(:,r),cfnParams.whatLoss,cfnParams.classLabels);
    end
end

%-------------------------------------------------------------------------------
%% Summarize quantities as the mean across repeats
%-------------------------------------------------------------------------------
if doSummarize
    accuracyMean = mean(accuracyMean,1);
    accuracyStd = mean(accuracyStd,1);
    inSampleAccStd = mean(inSampleAccStd,1);
end

if numFolds==0
    partitions = {};
end
cvOut = struct('partitions',{partitions},'transforms',{transforms},'yPredict',yPredict,...
                'numFeatures',numFeatures,'numImputed',numImputed);

end
