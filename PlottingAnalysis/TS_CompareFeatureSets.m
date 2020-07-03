function TS_CompareFeatureSets(whatData,whatClassifier,whatFeatureSets)
% TS_CompareFeatureSets Compares classification performance of feature sets
%
% Gives information about how different subsets of features behave on the data
% (length-dependent, location-dependent, spread-dependent features, and features
% that operate on the raw (rather than z-scored) time series)
%
% Runs a given classifier on the group labels assigned to the data, using
% different filters on the features.
%
% Provides a quick way of determining if there are location/spread/etc.
% differences between groups in a dataset
%
%---INPUTS:
% whatData: the dataset to analyze (input to TS_LoadData)
% whatClassifier: the classifier to apply to the different filters
% whatFeatureSets: custom set of feature-sets to compare against
%
%---USAGE:
% TS_CompareFeatureSets('norm','svm_linear');

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
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
if nargin < 2
    whatClassifier = 'svm_linear';
end
if nargin < 3
    whatFeatureSets = {'all','catch22','notLocationDependent','locationDependent',...
                        'notLengthDependent','lengthDependent',...
                        'notSpreadDependent','spreadDependent'};
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,dataFile] = TS_LoadData(whatData);
TS_GetFromData(whatData,'groupNames');

% Check that group labels have been assigned
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end

numFeatures = height(Operations);

%-------------------------------------------------------------------------------
% Define the feature sets by feature IDs
%-------------------------------------------------------------------------------
numFeatureSets = length(whatFeatureSets);
featureIDs = cell(numFeatureSets,1);
% Prep for pulling out IDs efficiently
opStruct = makeOpStruct();

for i = 1:numFeatureSets
    switch whatFeatureSets{i}
        case 'all'
            featureIDs{i} = Operations.ID;
        case 'notLengthDependent'
            [~,featureIDs{i}] = TS_GetIDs('lengthdep',opStruct,'ops','Keywords');
        case 'lengthDependent'
            featureIDs{i} = TS_GetIDs('lengthdep',opStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('lengthDependent',opStruct,'ops','Keywords');
        case 'notLocationDependent'
            [~,featureIDs{i}] = TS_GetIDs('locdep',opStruct,'ops','Keywords');
        case 'locationDependent'
            featureIDs{i} = TS_GetIDs('locdep',opStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('locationDependent',opStruct,'ops','Keywords');
        case 'notSpreadDependent'
            [~,featureIDs{i}] = TS_GetIDs('spreaddep',opStruct,'ops','Keywords');
        case 'spreadDependent'
            featureIDs{i} = TS_GetIDs('spreaddep',opStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('spreadDependent',opStruct,'ops','Keywords');
        case {'catch22','sarab16'}
            featureIDs{i} = GiveMeFeatureSet(whatFeatureSets{i},Operations);
        otherwise
            error('Unknown feature set: ''%s''',whatFeatureSets{i});
    end
end

numFeaturesIncluded = cellfun(@length,featureIDs);

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance
numClasses = max(TimeSeries.Group); % assumes group in form of integer class labels starting at 1
numFolds = HowManyFolds(TimeSeries.Group,numClasses);

fprintf(1,['Training and evaluating a %u-class %s classifier',...
                ' using %u-fold cross validation...\n'],...
                    numClasses,whatClassifier,numFolds);

accuracy = zeros(numFeatureSets,numFolds);
for i = 1:numFeatureSets
    filter = ismember(Operations.ID,featureIDs{i});
    [foldLosses,~,whatLoss] = GiveMeCfn(whatClassifier,TS_DataMat(:,filter),...
                    TimeSeries.Group,[],[],numClasses,[],[],true,numFolds,true);

    accuracy(i,:) = foldLosses;
    fprintf('Classified using the ''%s'' feature set (%u): (fold-average) %s = %.2f%%\n',...
                whatFeatureSets{i},numFeaturesIncluded(i),whatLoss,mean(accuracy(i,:)));
end


%-------------------------------------------------------------------------------
% Plot the result
dataCell = mat2cell(accuracy,ones(numFeatureSets,1),numFolds);
BF_JitteredParallelScatter(dataCell,true,true,true);
ax = gca();
ax.XTick = 1:numFeatureSets;
ax.XTickLabel = whatFeatureSets;
title(sprintf(['%u-class classification with different feature sets',...
                    ' using %u-fold cross validation'],...
                    numClasses,numFolds))
ax.TickLabelInterpreter = 'none';
ylabel(whatLoss)

%-------------------------------------------------------------------------------
function opStruct = makeOpStruct()
    % Generate a structure that contains the Operations table
    opStruct = struct();
    opStruct.TimeSeries = [];
    opStruct.TS_DataMat = [];
    opStruct.Operations = Operations;
end
%-------------------------------------------------------------------------------

end
