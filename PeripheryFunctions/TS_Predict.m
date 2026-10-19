function [tab,myAcc] = TS_Predict(timeSeriesData,timeSeriesNames,fileName_classifier,varargin)
% TS_Predict Predict classes for new time-series data using a pre-trained model
%
% This function uses a classifier saved by TS_Classify (and/or TS_TopFeatures)
% to predict the class labels of new time series.
% Only the features used by the classifier are computed, and they are transformed
% using the normalization parameters fitted on the training data (stored by
% TS_Normalize), so each time series is predicted as the trained model would
% predict it, independently of the other new time series (a single time series
% can be predicted).
%
%---USAGE:
% tab = TS_Predict(timeSeriesData,timeSeriesNames,fileName_classifier);
%
%---INPUTS:
% timeSeriesData, the new time-series data to predict: a cell array with one
%                 time series per cell (or a single time series as a vector)
% timeSeriesNames, a unique name for each time series (cell array of character
%                 vectors), or [] to name them ts_1, ts_2, ...
% fileName_classifier, the MAT-file previously saved using TS_Classify
%                   and/or TS_TopFeatures (cfnParams.classifierFilename)
%
% (OPTIONAL name-value pairs):
% classifierType, 'allFeatures' (the classifier from TS_Classify) or 'topFeature'
%                   (the single-feature classifier from TS_TopFeatures).
%                   Default: the one with the higher cross-validated accuracy if
%                   the file contains both.
% predictionFilename, a MAT-file to save the computed features (hctsa format) and
%                   the predictions to (default: '' -- not saved)
% overwrite [false], whether to overwrite an existing predictionFilename
% featureSet ['hctsa'], the feature set (input to TS_Init) that contains the
%                   classifier's features (matched by name)
% isParallel [false], how to parallelize the feature computation (input
%                   parallelHow to TS_Compute)
%
%---OUTPUTS:
% tab, a table with a row for each time series: Name, timeSeries (the data), and
%       predictedClass (categorical; undefined if a feature needed by the
%       classifier could not be computed for that time series)
% myAcc, the cross-validated performance of the classifier used (NaN if it was
%       not cross-validated)

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Oliver M. Cliff <oliver.m.cliff@gmail.com>,
% <http://www.olivercliff.github.io>
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
%% Check inputs
%-------------------------------------------------------------------------------
if isnumeric(timeSeriesData)
    timeSeriesData = {timeSeriesData}; % a single time series
end
timeSeriesData = timeSeriesData(:);
numTimeSeries = length(timeSeriesData);

if isempty(timeSeriesNames)
    timeSeriesNames = arrayfun(@(x)sprintf('ts_%u',x),(1:numTimeSeries)','UniformOutput',false);
end
timeSeriesNames = cellstr(timeSeriesNames);
timeSeriesNames = timeSeriesNames(:);
if length(timeSeriesNames)~=numTimeSeries
    error('%u time series but %u names',numTimeSeries,length(timeSeriesNames));
end
if length(unique(timeSeriesNames)) < numTimeSeries
    error('The names of the time series must be unique');
end

% Use an inputParser to control additional options as parameters:
inputP = inputParser;
addParameter(inputP,'classifierType','',@(x) ismember(x,{'','allFeatures','topFeature'}));
addParameter(inputP,'predictionFilename','',@ischar);
addParameter(inputP,'overwrite',false,@(x) islogical(x) || isnumeric(x));
addParameter(inputP,'featureSet','hctsa',@(x) ischar(x) || iscell(x));
addParameter(inputP,'isParallel',false,@(x) islogical(x) || ischar(x));
parse(inputP,varargin{:});
classifierType = inputP.Results.classifierType;
predictionFilename = inputP.Results.predictionFilename;
overwrite = inputP.Results.overwrite;
featureSet = inputP.Results.featureSet;
isParallel = inputP.Results.isParallel;
clear('inputP');

if ~isempty(predictionFilename) && exist(predictionFilename,'file') && ~overwrite
    error(['Prediction file ''%s'' already exists: set ''overwrite'' to true to replace it,' ...
                ' or set a different predictionFilename'],predictionFilename);
end

%-------------------------------------------------------------------------------
%% Load the classifier
%-------------------------------------------------------------------------------
fprintf('Loading pre-trained classifier from %s...', fileName_classifier);
classifierData = load(fileName_classifier); % jointClassifier and/or featureClassifier
fprintf('Done.\n');

hasJoint = isfield(classifierData,'jointClassifier') && ~isempty(classifierData.jointClassifier);
hasFeature = isfield(classifierData,'featureClassifier') && ~isempty(classifierData.featureClassifier) ...
                && isfield(classifierData.featureClassifier,'Mdl');
switch classifierType
    case 'allFeatures'
        useAllFeatures = true;
    case 'topFeature'
        useAllFeatures = false;
    otherwise
        % Select the model with the best cross-validated accuracy
        useAllFeatures = hasJoint && ~(hasFeature && ...
            classifierData.featureClassifier.CVAccuracy > classifierData.jointClassifier.CVAccuracy);
end
if useAllFeatures && ~hasJoint
    error('No joint classifier (from TS_Classify) in %s',fileName_classifier);
elseif ~useAllFeatures && ~hasFeature
    error('No single-feature classifier (from TS_TopFeatures) in %s',fileName_classifier);
end

if useAllFeatures
    % The joint classifier combines all features in the feature matrix
    myClassifier = classifierData.jointClassifier;
else
    % Otherwise, we use an individual classifier
    myClassifier = classifierData.featureClassifier;
end
myMdl = myClassifier.Mdl;
featureNames = cellstr(myClassifier.Operation.Name);
featureNames = featureNames(:); % the model's features, in the order of its columns
numFeatures = length(featureNames);
if isfield(myClassifier,'CVAccuracy')
    myAcc = myClassifier.CVAccuracy;
else
    myAcc = NaN;
end

if isa(myMdl,'classreg.learning.partition.ClassificationPartitionedModel')
    error(['The classifier in %s is a cross-validation object, not a trained classifier' ...
            ' (saved by an earlier version of TS_Classify): re-run TS_Classify to save' ...
            ' a usable classifier'],fileName_classifier);
end

if useAllFeatures
    fprintf('Using joint classifier (%u features) [CV acc=%.2f%%, in-sample acc=%.2f%% for %u classes].\n',...
            numFeatures,myAcc,myClassifier.Accuracy,length(myMdl.ClassNames));
else
    fprintf('Using classifier for feature "%s" [CV acc=%.2f%%, in-sample acc=%.2f%% for %u classes].\n',...
            featureNames{1},myAcc,myClassifier.Accuracy,length(myMdl.ClassNames));
end

%-------------------------------------------------------------------------------
%% Get the normalizing transformation of the training data
%-------------------------------------------------------------------------------
if isfield(myClassifier,'normalizationInfo')
    normInfo = myClassifier.normalizationInfo;
else
    normInfo = [];
end
doNormalize = ~isempty(normInfo) && ~ismember(normInfo.normFunction,{'none','nothing'});
if doNormalize
    if ~isfield(normInfo,'featureParams')
        error(['The classifier was trained on data normalized (''%s'') by a version of' ...
                ' TS_Normalize that did not store the normalization parameters, so new data' ...
                ' cannot be transformed in the same way. Re-run TS_Normalize and then' ...
                ' TS_Classify (or TS_TopFeatures) to save a usable classifier.'],normInfo.normFunction);
    end
    [hasParams,ind] = ismember(featureNames,normInfo.featureParams.Name);
    if ~all(hasParams)
        error('No normalization parameters stored for %u features of the classifier (e.g., %s)',...
                    sum(~hasParams),featureNames{find(~hasParams,1)});
    end
    normParams = struct('normMethod',normInfo.normFunction,...
                        'featureParams',normInfo.featureParams(ind,:));
end

%-------------------------------------------------------------------------------
%% Compute the features needed by the classifier
%-------------------------------------------------------------------------------
% Work in a new temporary directory (removed on exit, including on error):
workDir = tempname;
mkdir(workDir);
cleanupObj = onCleanup(@() rmdir(workDir,'s'));
inputFile = fullfile(workDir,'INP_ts.mat');
hctsaFile = fullfile(workDir,'HCTSA.mat');

labels = timeSeriesNames;
keywords = repmat({''},numTimeSeries,1);
save(inputFile,'timeSeriesData','labels','keywords');
TS_Init(inputFile,featureSet,false,hctsaFile);

% Find the classifier's features in the feature set (by name):
Operations = TS_GetFromData(hctsaFile,'Operations');
[isFound,opInd] = ismember(featureNames,Operations.Name);
if ~all(isFound)
    error(['%u features of the classifier are not in the feature set (e.g., %s):' ...
            ' set featureSet to the feature set the classifier was trained on'],...
            sum(~isFound),featureNames{find(~isFound,1)});
end
opIDs = Operations.ID(opInd);

% Compute just these features:
TS_Compute(isParallel,[],opIDs,'missing',hctsaFile,'minimal');

% Load the computed features (bad values as NaN), in the order of the model's columns:
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(hctsaFile);
[~,rowInd] = ismember(timeSeriesNames,TimeSeries.Name);
[~,colInd] = ismember(opIDs,Operations.ID);
featureMatrix = TS_DataMat(rowInd,colInd);

%-------------------------------------------------------------------------------
%% Transform the features as the training data, and predict
%-------------------------------------------------------------------------------
if doNormalize
    featureMatrix = BF_NormalizeMatrix(featureMatrix,normParams);
end

% Only predict time series with all the features needed:
isGood = all(isfinite(featureMatrix),2);
if ~all(isGood)
    warning(['%u time series have special values for features needed by the classifier,' ...
                ' so are not predicted: %s'],sum(~isGood),BF_cat(timeSeriesNames(~isGood),','));
end
classNames = cellstr(myMdl.ClassNames);
predictedClass = categorical(repmat({''},numTimeSeries,1),classNames); % (undefined)
if any(isGood)
    predictedClass(isGood) = cellstr(predict(myMdl,featureMatrix(isGood,:)));
end

tab = table(timeSeriesNames,timeSeriesData,predictedClass,'VariableNames',...
                {'Name','timeSeries','predictedClass'});

%-------------------------------------------------------------------------------
%% Save to file
%-------------------------------------------------------------------------------
if ~isempty(predictionFilename)
    % The computed features (hctsa format), the transformed features input to
    % the classifier (in the order of featureNames), and the predictions:
    copyfile(hctsaFile,predictionFilename);
    predictionTable = tab;
    save(predictionFilename,'predictionTable','featureMatrix','featureNames','-append');
    fprintf('Saved computed features and predictions to %s\n',predictionFilename);
end

end
