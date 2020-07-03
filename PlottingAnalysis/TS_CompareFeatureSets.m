function TS_CompareFeatureSets(whatData,whatFeatureSets,cfnParams)
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
% differences between groups in a dataset.
%
%---INPUTS:
% whatData: the dataset to analyze (input to TS_LoadData)
% whatFeatureSets: custom set of feature-sets to compare against
% cfnParams: custom classification parameters for running classificaiton algorithms
%
%---USAGE:
% TS_CompareFeatureSets('norm');

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
    whatFeatureSets = {'all','catch22','notLocationDependent','locationDependent',...
                        'notLengthDependent','lengthDependent',...
                        'notSpreadDependent','spreadDependent'};
end
% <set default classificaiton parameters after inspecting TimeSeries labeling below>

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,dataFile] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
dataStruct = makeDataStruct();
numFeatures = height(Operations);

TellMeAboutLabeling(dataStruct);

%-------------------------------------------------------------------------------
% Set classification parameters
%-------------------------------------------------------------------------------
if nargin < 3
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end

%-------------------------------------------------------------------------------
% Define the feature sets by feature IDs
%-------------------------------------------------------------------------------
numFeatureSets = length(whatFeatureSets);
featureIDs = cell(numFeatureSets,1);
theColors = cell(numFeatureSets,1);
% Prep for pulling out IDs efficiently

for i = 1:numFeatureSets
    switch whatFeatureSets{i}
        case 'all'
            featureIDs{i} = Operations.ID;
            theColors{i} = [203,90,76]/255;
        case 'notLengthDependent'
            [~,featureIDs{i}] = TS_GetIDs('lengthdep',dataStruct,'ops','Keywords');
            theColors{i} = brighten([193,92,165]/255,+0.5);
        case 'lengthDependent'
            featureIDs{i} = TS_GetIDs('lengthdep',dataStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('lengthDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([193,92,165]/255,-0.5);
        case 'notLocationDependent'
            [~,featureIDs{i}] = TS_GetIDs('locdep',dataStruct,'ops','Keywords');
            theColors{i} = brighten([180,148,62]/255,+0.5);
        case 'locationDependent'
            featureIDs{i} = TS_GetIDs('locdep',dataStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('locationDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([180,148,62]/255,-0.5);
        case 'notSpreadDependent'
            [~,featureIDs{i}] = TS_GetIDs('spreaddep',dataStruct,'ops','Keywords');
            theColors{i} = brighten([114,124,206]/255,+0.5);
        case 'spreadDependent'
            featureIDs{i} = TS_GetIDs('spreaddep',dataStruct,'ops','Keywords');
            % featureIDs{i} = TS_GetIDs('spreadDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([114,124,206]/255,-0.5);
        case {'catch22','sarab16'}
            featureIDs{i} = GiveMeFeatureSet(whatFeatureSets{i},Operations);
            theColors{i} = [96,168,98]/255;
        otherwise
            error('Unknown feature set: ''%s''',whatFeatureSets{i});
    end
end

numFeaturesIncluded = cellfun(@length,featureIDs);

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance
% This needs to be on for the below to work:
cfnParams.computePerFold = true;
TellMeAboutClassification(cfnParams);

accuracy = zeros(numFeatureSets,cfnParams.numFolds*cfnParams.numRepeats);
for i = 1:numFeatureSets
    filter = ismember(Operations.ID,featureIDs{i});
    for j = 1:cfnParams.numRepeats
        [foldLosses,~,whatLoss] = GiveMeCfn(cfnParams.whatClassifier,TS_DataMat(:,filter),...
                    TimeSeries.Group,[],[],cfnParams.numClasses,[],[],cfnParams.doReweight,cfnParams.numFolds,true);
        accuracy(i,1+(j-1)*cfnParams.numFolds:j*cfnParams.numFolds) = foldLosses;
    end
    fprintf(['Classified using the ''%s'' set (%u features): (%u fold-average, ',...
                        '%u repeats) average %s = %.2f%%\n'],...
            whatFeatureSets{i},numFeaturesIncluded(i),cfnParams.numFolds,cfnParams.numRepeats,...
            whatLoss,mean(accuracy(i,:)));
end


%-------------------------------------------------------------------------------
% Plot the result
dataCell = mat2cell(accuracy,ones(numFeatureSets,1),size(accuracy,2));
extraParams = struct();
extraParams.theColors = theColors;
BF_JitteredParallelScatter(dataCell,true,true,true,extraParams);
ax = gca();
ax.XTick = 1:numFeatureSets;
ax.XTickLabel = whatFeatureSets;
ax.XTickLabelRotation = 45;
title(sprintf(['%u-class classification with different feature sets',...
                    ' using %u-fold cross validation'],...
                    cfnParams.numClasses,cfnParams.numFolds))
ax.TickLabelInterpreter = 'none';
ylabel(whatLoss)

%-------------------------------------------------------------------------------
function dataStruct = makeDataStruct()
    % Generate a structure for the dataset
    dataStruct = struct();
    dataStruct.TimeSeries = TimeSeries;
    dataStruct.TS_DataMat = TS_DataMat;
    dataStruct.Operations = Operations;
    dataStruct.groupNames = TS_GetFromData(whatData,'groupNames');
end
%-------------------------------------------------------------------------------

end
