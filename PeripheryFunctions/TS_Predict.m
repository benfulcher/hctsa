function [tab,myAcc] = TS_Predict(timeSeriesData,labels,fileName_classifier,varargin)
% TS_predict Predict classes using new time-series data
%
% This function uses a previously learnt classifier to predict group labels
% of input time series.
%
%---USAGE:
% TS_predict(timeSeriesData,labels,fileName_classifier);
%
%---INPUTS:
% timeSeriesData, the new time-series data to predict (cell array)
% labels, the labels for each time series
% fileName_classifier, the MAT-file previously saved using TS_classify
% and/or TS_TopFeatures
%
% (OPTIONAL):
% classifierType, 'topFeature', or 'allFeatures'
% predictionFilename, the output filename containing predicted keywords,
% and associated input time series/labels
%
%---OUTPUTS:
% Table of timeSeries, labels, and their predicted keywords (classes)

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

% Use an inputParser to control additional plotting options as parameters:
inputP = inputParser;
% classifierType:
default_classifierType = '';
check_classifierType = @(x) ischar(x);
addParameter(inputP,'classifierType',default_classifierType,check_classifierType);

default_predictionFilename = '';
check_predictionFilename = @(x) ischar(x);
addParameter(inputP,'predictionFilename',default_predictionFilename,check_predictionFilename);

default_isParallel = '';
check_isParallel = @(x) islogical(x);
addParameter(inputP,'isParallel',default_isParallel,check_isParallel);

% Parse input arguments:
parse(inputP,varargin{:});
classifierType = inputP.Results.classifierType;
predictionFilename = inputP.Results.predictionFilename;
isParallel = inputP.Results.isParallel;
clear('inputP');

%#ok<*NASGU>

fprintf('Loading %s...', fileName_classifier);
load(fileName_classifier);
fprintf('Done.\n');

classes = reshape(classes,length(classes),1);

%-------------------------------------------------------------------------------
% Choose classifier:
%-------------------------------------------------------------------------------
useAllFeatures = true;
if ~isempty(classifierType)
    % If we've explicitly asked for a classifier type
    switch classifierType
        case 'allFeatures'
            useAllFeatures = true;
        case 'topFeature'
            useAllFeatures = false;
    end
else
    % Otherwise, select the model with the best cross-validated accuracy
    if exist('featureClassifier','var')
      if exist('jointClassifier','var')
        if any(featureClassifier.CVAccuracy > jointClassifier.CVAccuracy)
          useAllFeatures = false;
        end
      else
        useAllFeatures = false;
      end
    elseif ~exist('jointClassifier','var')
      useAllFeatures = false;
    end 
end

if useAllFeatures
  myClassifier = jointClassifier;
  myAcc = myClassifier.Accuracy;
  fprintf('Using joint classifier (%i features) [acc=%.2f%% for %i classes].\n',...
                length(jointClassifier.Operation.ID),...
                myAcc, length(myClassifier.Mdl.ClassNames));
else
  myClassifier = featureClassifier;
  myAcc = myClassifier.Accuracy;
  
  fprintf('Using classifier for feature "%s" (%i) [acc=%.2f%% for %i classes].\n',...
                myClassifier.Operation.Name, myClassifier.Operation.ID,...
                myAcc, length(myClassifier.Mdl.ClassNames));
end

myOps = myClassifier.Operation.ID;
myMdl = myClassifier.Mdl;

removeFile = false;

toRecompute = true;
if exist(predictionFilename,'file')
    out = input(sprintf('Warning: %s already exists -- override? [yn] ',...
                    predictionFilename),'s');
    if out ~= 'y'
        toRecompute = false;
    end
end
    
if toRecompute
    if isempty(predictionFilename)
        removeFile = true;
        predictionFilename = 'tmp';
        i = 1;
        while exist([predictionFilename '.mat'],'file')
            predictionFilename = sprintf('tmp-%i',i);
            i = i + 1;
        end
        predictionFilename = [predictionFilename '.mat'];
    end

    keywords = cell(length(timeSeriesData),1);

    [filePath,name] = fileparts(predictionFilename);

    if ~isempty(filePath)
        tsFilename = [filePath '/' name '_T.mat'];
    else
        tsFilename = [name '_T.mat'];
    end

    if exist(tsFilename,'file')
        [~,name] = fileparts(tsFilename);
        out = input(sprintf('Warning: %s already exists -- override? [yn] ',...
                    name),'s');
        if out ~= 'y'
            return
        end
    end

    save(tsFilename,'timeSeriesData','labels','keywords','-v7.3');
    TS_Init(tsFilename,'','',0,predictionFilename);
    delete(tsFilename);

    TS_Compute(isParallel,[],myOps,'',predictionFilename,0);
end

if isfield(myClassifier,'normalizationInfo')
    TS_Normalize(myClassifier.normalizationInfo.normFunction,...
                    [0 1],...
                    predictionFilename);
                
    if exist([predictionFilename(1:end-4) '_N.mat'],'file')
        delete(predictionFilename);
        predictionFilename = [predictionFilename(1:end-4) '_N.mat'];
    end
end

whatData = load(predictionFilename);

[TS_DataMat,~,Operations] = TS_LoadData(whatData);

[myOpId,~] = find(Operations.ID == myOps');

predictGroups = predict(myMdl,TS_DataMat(:,myOpId));
predictKeywords = classes(predictGroups); % Predict labels from new time-series data

tab = table(labels,timeSeriesData,predictKeywords,predictGroups,'VariableNames',{'labels','timeSeries','predictKeywords','predictGroups'});

if removeFile
    delete(predictionFilename);
else
    save(predictionFilename,'predictGroups','predictKeywords','-append');
end