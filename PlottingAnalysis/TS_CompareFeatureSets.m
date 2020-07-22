function TS_CompareFeatureSets(whatData,cfnParams,whatFeatureSets)
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
% cfnParams: custom classification parameters for running classificaiton algorithms
% whatFeatureSets: custom set of feature-sets to compare against
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
if nargin < 3 || isempty(whatFeatureSets)
    whatFeatureSets = {'all','catch22','notLocationDependent','locationDependent',...
                        'notLengthDependent','lengthDependent',...
                        'notSpreadDependent','spreadDependent'};
end
% <set default classificaiton parameters after inspecting TimeSeries labeling below>

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,dataFile] = TS_LoadData(whatData);

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
[groupLabels,classLabels,groupLabelsInteger,numGroups] = ExtractGroupLabels(TimeSeries);
TellMeAboutLabeling(TimeSeries);

% Set classification parameters if needed:
if nargin < 2 || isempty(cfnParams)
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end

% Set up data for Get_IDs
dataStruct = makeDataStruct();

%-------------------------------------------------------------------------------
% Define the feature sets as sets of IDs
%-------------------------------------------------------------------------------
numFeatureSets = length(whatFeatureSets);
featureIDs = cell(numFeatureSets,1);
theColors = cell(numFeatureSets,1);
featureSetNames = cell(numFeatureSets,1);
% Prep for pulling out IDs efficiently

for i = 1:numFeatureSets
    switch whatFeatureSets{i}
        case 'all'
            featureIDs{i} = Operations.ID;
            featureSetNames{i} = sprintf('hctsa (%u)',height(Operations));
            theColors{i} = [233,129,126]/255;
        case {'catch22','sarab16'}
            featureIDs{i} = GiveMeFeatureSet(whatFeatureSets{i},Operations);
            featureSetNames{i} = sprintf('%s (%u)',whatFeatureSets{i},length(featureIDs{i}));
            theColors{i} = [151,205,104]/255;
        case 'notLengthDependent'
            [~,featureIDs{i}] = TS_GetIDs('lengthDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('hctsa without length-dependent (%u)',length(featureIDs{i}));
            theColors{i} = brighten([187,149,219]/255,+0.3);
        case 'lengthDependent'
            featureIDs{i} = TS_GetIDs('lengthDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('Length-dependent (%u)',length(featureIDs{i}));
            % featureIDs{i} = TS_GetIDs('lengthDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([187,149,219]/255,-0.3);
        case 'notLocationDependent'
            [~,featureIDs{i}] = TS_GetIDs('locationDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('hctsa without location-dependent (%u)',length(featureIDs{i}));
            theColors{i} = brighten([214,175,90]/255,+0.3);
        case 'locationDependent'
            featureIDs{i} = TS_GetIDs('locationDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('Location-dependent (%u)',length(featureIDs{i}));
            % featureIDs{i} = TS_GetIDs('locationDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([214,175,90]/255,-0.3);
        case 'notSpreadDependent'
            [~,featureIDs{i}] = TS_GetIDs('spreadDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('hctsa without spread-dependent (%u)',length(featureIDs{i}));
            theColors{i} = brighten([111,204,180]/255,+0.3);
        case 'spreadDependent'
            featureIDs{i} = TS_GetIDs('spreadDependent',dataStruct,'ops','Keywords');
            featureSetNames{i} = sprintf('Spread-dependent (%u)',length(featureIDs{i}));
            % featureIDs{i} = TS_GetIDs('spreadDependent',dataStruct,'ops','Keywords');
            theColors{i} = brighten([111,204,180]/255,-0.3);
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
        foldLosses = GiveMeCfn(TS_DataMat(:,filter),TimeSeries.Group,[],[],cfnParams);
        accuracy(i,1+(j-1)*cfnParams.numFolds:j*cfnParams.numFolds) = foldLosses;
    end
    fprintf(['Classified using the ''%s'' set (%u features): (%u-fold average, ',...
                        '%u repeats) average %s = %.2f%%\n'],...
            whatFeatureSets{i},numFeaturesIncluded(i),cfnParams.numFolds,cfnParams.numRepeats,...
            cfnParams.whatLoss,mean(accuracy(i,:)));
end


%-------------------------------------------------------------------------------
% Plot the result
dataCell = mat2cell(accuracy,ones(numFeatureSets,1),size(accuracy,2));
extraParams = struct();
extraParams.theColors = theColors;
f = figure('color','w');
% Add clear horizontal comparison to performance of all features
if ismember('all',whatFeatureSets)
    allIndex = find(strcmp(whatFeatureSets,'all'));
    plot([1,numFeatureSets],mean(dataCell{allIndex}),'--','color',theColors{allIndex})
end
BF_JitteredParallelScatter(dataCell,true,true,false,extraParams);
ax = gca();
ax.XTick = 1:numFeatureSets;
ax.XTickLabel = whatFeatureSets;
ax.XTickLabelRotation = 30;
title(sprintf(['%u-class classification with different feature sets',...
                    ' using %u-fold cross validation'],...
                    cfnParams.numClasses,cfnParams.numFolds))
ax.TickLabelInterpreter = 'none';
ylabel(sprintf('%s (%s)',cfnParams.whatLoss,cfnParams.whatLossUnits))
f.Position(3:4) = [732   423];


%-------------------------------------------------------------------------------
function dataStruct = makeDataStruct()
    % Generate a structure for the dataset
    dataStruct = struct();
    dataStruct.TimeSeries = TimeSeries;
    dataStruct.TS_DataMat = TS_DataMat;
    dataStruct.Operations = Operations;
end
%-------------------------------------------------------------------------------

end
