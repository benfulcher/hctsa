function TS_compareFilters(whatData,whatClassifier)
% TS_compareFilters  Gives information about how different subsets of features
%                   behave on the data (length-dependent, location-dependent,
%                   spread-dependent features, and features operating on
%                   the raw, rather than z-scored, time series)
%
% Runs a given classifier on the group labels assigned to the data, using
% different filters on the features.
% Provides a quick way of determining if there are location/spread/etc.
% differences between groups in a dataset
%
%---INPUTS:
% whatData: the dataset to analyse (input to TS_LoadData)
% whatClassifier: the classifier to apply to the different filters
%
%---USAGE:
% TS_compareFilters('norm','svm_linear');

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,dataFile] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~isfield(TimeSeries,'Group')
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
load(dataFile,'groupNames');
timeSeriesGroup = [TimeSeries.Group]'; % Use group form (column vector)
numClasses = max(timeSeriesGroup); % assuming group in form of integer class labels starting at 1
numFeatures = length(Operations);
numFolds = howManyFolds(timeSeriesGroup,numClasses);

%-------------------------------------------------------------------------------
% Get the filters
%-------------------------------------------------------------------------------
dataStruct.TimeSeries = [];
dataStruct.TS_DataMat = [];
dataStruct.Operations = Operations;
[ID_lengthDep,ID_notlengthDep] = TS_getIDs('lengthdep',dataStruct,'ops');
[ID_locDep,ID_notlocDep] = TS_getIDs('locdep',dataStruct,'ops');
[ID_spreadDep,ID_notspreadDep] = TS_getIDs('spreaddep',dataStruct,'ops');
[ID_raw,ID_notraw] = TS_getIDs('raw',dataStruct,'ops');

%-------------------------------------------------------------------------------
% Compare each
%-------------------------------------------------------------------------------
accuracy = zeros(6,numFolds);
numFeaturesIncluded = zeros(6,1);
names = cell(6,1);
% (i) all features, (ii) no length, (iii) no spread, (iv) no location, (v) no length/spread/location

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance

fprintf(1,['Training and evaluating a %u-class %s classifier using %u-fold cross validation...\n'],...
                        numClasses,whatClassifier,numFolds);

for i = 1:6
    switch i
    case 1
        filter = true(length(Operations),1);
        names{1} = 'all';
    case 2
        filter = ismember([Operations.ID],ID_notlengthDep);
        names{2} = 'no length';
    case 3
        filter = ismember([Operations.ID],ID_notlocDep);
        names{3} = 'no location';
    case 4
        filter = ismember([Operations.ID],ID_notspreadDep);
        names{4} = 'no spread';
    case 5
        filter = ismember([Operations.ID],ID_notraw);
        names{5} = 'no raw';
    case 6
        filter = ismember([Operations.ID],intersect(intersect(ID_notlengthDep, ID_notlocDep), ID_notspreadDep));
        names{6} = 'no length, location, spread';
    end
    [foldLosses,~,whatLoss] = GiveMeCfn(whatClassifier,TS_DataMat(:,filter),...
                            timeSeriesGroup,[],[],numClasses,[],[],1,numFolds);

    accuracy(i,:) = foldLosses;
    numFeaturesIncluded(i) = sum(filter);
    names{i} = sprintf('%s (%u)',names{i},numFeaturesIncluded(i));
end

%-------------------------------------------------------------------------------
% Plot the result
f = figure('color','w'); ax = gca;
errorbar(1:6,mean(accuracy,2),std(accuracy,[],2))
ax.XTick = 1:6;
ax.XTickLabel = names;
title(sprintf('%u-class classification with different feature sets using %u-fold cross validation',...
                        numClasses,numFolds))
ax.TickLabelInterpreter = 'none';
ylabel(whatLoss)


end
