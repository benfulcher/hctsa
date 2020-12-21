function params = GiveMeDefaultClassificationParams(TimeSeries,numClasses,beVocal)
% Set default parameters describing the feature-based time-series classification
%

%-------------------------------------------------------------------------------
% If you use this code for your research, please cite these papers:
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
%-------------------------------------------------------------------------------

% Get TimeSeries labeling information from HCTSA.mat by default
if nargin < 1 || isempty(TimeSeries)
    [~,TimeSeries] = TS_LoadData('HCTSA.mat');
    warning('DEFAULT: Retrieving time-series labeling information from HCTSA.mat')
end
if ~istable(TimeSeries)
    % Actually specified TimeSeries as a filename or structure:
    TimeSeries = TS_GetFromData(TimeSeries,'TimeSeries');
end
if nargin < 3
    beVocal = true;
end
%-------------------------------------------------------------------------------
% Set the classifier:
params.whatClassifier = 'svm_linear'; % ('svm_linear', 'knn', 'linear', 'fast_linear')

% Number of repeats of cross-validation (reduce variance due to 'lucky splits')
params.numRepeats = 1;

% Check group labeling:
if ~ismember('Group',TimeSeries.Properties.VariableNames) || all(isundefined(TimeSeries.Group))
    warning('Group labels not assigned to time series. Cannot perform classification.');
    params = struct();
    return
end

% Number of classes to classify
% (Assume every class is represented in the data):
params.classLabels = categories(TimeSeries.Group);
if nargin < 2 || isempty(numClasses);
    numClasses = length(params.classLabels);
end
params.numClasses = numClasses;

% Get numbers in each class:
classNumbers = arrayfun(@(x)sum(TimeSeries.Group==x),params.classLabels);
isBalanced = all(classNumbers==classNumbers(1));

% Set the performance metric, and balancing settings based on class balance statistics
if isBalanced
    params.doReweight = false;
    params.whatLoss = 'accuracy';
    params.whatLossUnits = '%';
else
    if beVocal
        fprintf(1,'Unbalanced classes: using a balanced accuracy measure (& using reweighting)...\n');
    end
    params.doReweight = true;
    params.whatLoss = 'balancedAccuracy';
    params.whatLossUnits = '%';
end

% Cross validation: number of folds (set to 0 for no CV)
params.numFolds = HowManyFolds(TimeSeries.Group,numClasses);

% Whether to output information about each fold, or average over folds
params.computePerFold = false;

% .mat file to save the classifier to (not saved if empty).
params.classifierFilename = ''; % (don't save classifier information to file)

% Set text to name the classifier
params = UpdateClassifierText(params);

% Restrict to a reduced set of features?
% '' (all), 'catch22', 'catchaMouse16', 'noLengthLocationSpread'
% (cf. TS_GiveMeFeatureSet)
params.reducedFeatureSet = '';

end
