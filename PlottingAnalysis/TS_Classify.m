function [meanAcc,nullStats,jointClassifier] = TS_Classify(whatData,cfnParams,numNulls,varargin)
% TS_Classify   Classify class labels assigned to the data using all features
%
% Requires the time-series data to have been labeled using TS_LabelGroups, which
% stores in the Group column of the TimeSeries table.
%
%---USAGE:
% TS_Classify();
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData)
% cfnParams, parameters for the classification algorithm (cf. GiveMeDefaultClassificationParams)
% numNulls (integer), number of nulls to compute (0 for no null comparison)
%
% (OPTIONAL):
% doParallel [false], whether to speed up the null computation using a parfor loop
% doPlot [true], whether to plot outputs (including confusion matrix)
% seedReset, input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)
% overwrite [false], whether to overwrite an existing cfnParams.classifierFilename
% normalizeWithinFolds [false], whether to fit the feature filter and normalization
%               within each training fold (see below): true (using the settings
%               recorded by TS_Normalize in the normalized data, or the defaults
%               of TS_Normalize for data that were not normalized), or a
%               normalization specification structure (cf. BF_FitTransform:
%               fields normFunction, filterOptions, classVarFilter, testNaN).
% rawData, the unnormalized data used when normalizeWithinFolds is set (input
%               to TS_LoadData). Default: the file TS_Normalize was given (the
%               name of the normalized file without '_N'), or whatData itself if
%               it is not normalized.
%
%---NORMALIZATION AND CROSS-VALIDATION:
% By default, the classifier is cross-validated on the data as given. For data
% normalized by TS_Normalize, this means that the feature filter and the
% normalization of each feature (e.g., its median and IQR) were fitted on all
% time series, including those in each test fold. No labels are used (unless
% TS_Normalize was run with classVarFilter), but the test folds are not
% entirely unseen, so cross-validated accuracies can be slightly optimistic.
% To avoid this, set normalizeWithinFolds: the unnormalized data are loaded
% (keeping the time series of whatData), and for each fold the feature filter
% and normalization are fitted on the training fold only (BF_FitTransform) and
% applied to the training and test folds (BF_ApplyTransform), before the
% classifier is trained on the training fold and evaluated on the test fold
% (BF_CrossValidate). Missing values of kept features in a test fold are
% replaced by the median of the normalized training-fold values (unless the
% specification sets testNaN = 'keep'); the number replaced is reported.
% Model-based nulls use the same partitions as the real labels, with a
% transformation fitted once per fold and reused across permutations (refitted
% for each permutation if classVarFilter is set, as it depends on the labels);
% label permutations are drawn with every class in every training and test fold.
%
%---SAVING A CLASSIFIER:
% If cfnParams.classifierFilename is set, a classifier is trained on all the
% data (no cross-validation) and saved to that file (replacing it if overwrite
% is true; an existing file gives an error otherwise), for use by TS_Predict.
% The file contains jointClassifier, classes (the class labels), and cfnParams.
%
%---OUTPUTS:
% (to commandline) output on classification performance on the dataset using the feature set.
% - meanAcc, the performance metric across repeats of cross-validation
% - nullStats, the performance metric across randomizations of the data labels
% - jointClassifier, the saved all-features classifier ([] if none is saved),
%       a structure with fields:
%       Mdl: the classification model, trained on all the data
%       Operation: the features used (ID and Name), in the order of the model's columns
%       classLabels: the class labels
%       Accuracy: in-sample performance of Mdl (cfnParams.whatLoss)
%       CVAccuracy: cross-validated performance (mean across folds and repeats;
%                   NaN if cfnParams.numFolds = 0)
%       CVAccuracyFolds: performance in each test fold (numFolds x numRepeats; empty
%                   unless cfnParams.computePerFold and cfnParams.numFolds > 0)
%       whatLoss, whatLossUnits, whatClassifier, numFolds, numRepeats: from cfnParams
%       normalizationInfo: how the data were normalized (from TS_Normalize), with
%                   featureParams restricted to the features of the model, in
%                   the same order ([] for data that were not normalized). With
%                   normalizeWithinFolds, the transformation is fitted on all
%                   the time series (BF_FitTransform), and the model uses the
%                   features it keeps.

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
%% Check Inputs
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
% <set default classification parameters after inspecting TimeSeries labeling below>
if nargin < 3
    numNulls = 0;
end

% Use an inputParser to parse additional options as parameters:
inputP = inputParser;

% simpleNull: use simple label shuffles as nulls
default_simpleNull = true;
check_simpleNull = @(x) islogical(x);
addParameter(inputP,'simpleNull',default_simpleNull,check_simpleNull);

% doParallel: use parfor to speed up null computation
default_doParallel = false;
check_doParallel = @(x) islogical(x);
addParameter(inputP,'doParallel',default_doParallel,check_doParallel);

% seedReset:
default_seedReset = 'default';
check_seedReset = @(x) ischar(x);
addParameter(inputP,'seedReset',default_seedReset,check_seedReset);

% doPlot:
default_doPlot = true;
check_doPlot = @(x) islogical(x);
addParameter(inputP,'doPlot',default_doPlot,check_doPlot);

% overwrite: whether to overwrite an existing classifier file
default_overwrite = false;
check_overwrite = @(x) islogical(x) || isnumeric(x);
addParameter(inputP,'overwrite',default_overwrite,check_overwrite);

% normalizeWithinFolds: fit the feature filter and normalization within each training fold
default_normalizeWithinFolds = false;
check_normalizeWithinFolds = @(x) islogical(x) || isstruct(x) || (isnumeric(x) && isscalar(x));
addParameter(inputP,'normalizeWithinFolds',default_normalizeWithinFolds,check_normalizeWithinFolds);

% rawData: the unnormalized data (used with normalizeWithinFolds)
default_rawData = '';
check_rawData = @(x) ischar(x) || isstring(x) || isstruct(x);
addParameter(inputP,'rawData',default_rawData,check_rawData);

% Parse input arguments:
parse(inputP,varargin{:});
simpleNull = inputP.Results.simpleNull;
doParallel = inputP.Results.doParallel;
seedReset = inputP.Results.seedReset;
doPlot = inputP.Results.doPlot;
overwrite = inputP.Results.overwrite;
normalizeWithinFolds = inputP.Results.normalizeWithinFolds;
rawData = inputP.Results.rawData;
clear('inputP');
doWithinFolds = isstruct(normalizeWithinFolds) || logical(normalizeWithinFolds);

%-------------------------------------------------------------------------------
%% Load in data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);

% Use the unnormalized data to normalize within folds:
if doWithinFolds
    [TS_DataMat,TimeSeries,Operations,normSpec] = LoadRawData(whatData,whatDataFile,...
                                    TimeSeries,rawData,normalizeWithinFolds);
end
[~,classLabels] = TS_ExtractGroupLabels(TimeSeries);

% Give basic info about the represented classes:
TellMeAboutLabeling(TimeSeries);

% Settings for the classification model:
if nargin < 2 || (isstruct(cfnParams) && isempty(fields(cfnParams))) || isempty(cfnParams)
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end
TellMeAboutClassification(cfnParams);

% Check whether the classifier can be saved before doing any work:
doSaveClassifier = isfield(cfnParams,'classifierFilename') && ~isempty(cfnParams.classifierFilename);
if doSaveClassifier && ~overwrite && exist(cfnParams.classifierFilename,'file')
    error(['Classifier file ''%s'' already exists: set ''overwrite'' to true to replace it,' ...
                ' or set a different cfnParams.classifierFilename'],cfnParams.classifierFilename);
end

% Filter down a reduced feature set if required:
[TS_DataMat,Operations] = FilterFeatures(TS_DataMat,Operations,cfnParams);
numFeatures = height(Operations);

%-------------------------------------------------------------------------------
%% Fit the model
%-------------------------------------------------------------------------------
% Reset the random seed for CV-reproducibility
BF_ResetSeed(seedReset);
fprintf(1,'Resetting random seed for reproducibility\n');

% Fit the classification model to the dataset and evaluate performance:
if doWithinFolds
    [meanAcc,stdAcc,inSampleAccStd,allFoldData,cvOut] = BF_CrossValidate(TS_DataMat,...
                                        TimeSeries.Group,cfnParams,normSpec,...
                                        'Operations',Operations,'doSummarize',false);
    featureText = sprintf('%u-%u features (filtered and normalized within each training fold, from %u)',...
                        min(cvOut.numFeatures(:)),max(cvOut.numFeatures(:)),numFeatures);
    fprintf(1,'%s normalization fitted within each training fold: %s.\n',normSpec.normFunction,featureText);
    if any(cvOut.numImputed(:) > 0)
        fprintf(1,['%u missing values of kept features in test folds (%.2g%% of test-fold values)' ...
                    ' were replaced by the training-fold median.\n'],sum(cvOut.numImputed(:)),...
                    100*sum(cvOut.numImputed(:))/sum(cvOut.numFeatures.*TestFoldSizes(cvOut,height(TimeSeries)),'all'));
    end
else
    [meanAcc,stdAcc,inSampleAccStd,allFoldData,CVMdl] = ComputeCVAccuracies(TS_DataMat,...
                                        TimeSeries.Group,cfnParams,false);
    featureText = sprintf('%u features',numFeatures);
end
% Whether results are available for each fold:
perFold = cfnParams.computePerFold && (cfnParams.numFolds > 0);

%-------------------------------------------------------------------------------
%% Display results to command line
%-------------------------------------------------------------------------------
if cfnParams.numFolds==0
    assert(cfnParams.numRepeats==1)
    fprintf(1,'In-sample %s using %s = %.3f%s\n',cfnParams.whatLoss,cfnParams.whatClassifier,...
                                            meanAcc,cfnParams.whatLossUnits);
else
    if cfnParams.numRepeats==1
        if cfnParams.computePerFold
            whatStat = sprintf('Mean +/- std (across %u folds)',cfnParams.numFolds);
            accuracyBit = sprintf('%.3f +/- %.3f%s',meanAcc,stdAcc,cfnParams.whatLossUnits);
        else
            whatStat = sprintf('Mean (across %u folds)',cfnParams.numFolds);
            accuracyBit = sprintf('%.3f%s',meanAcc,cfnParams.whatLossUnits);
        end
    else
        if cfnParams.computePerFold
            whatStat = sprintf('Mean +/- std (across %u folds, means across %u repeats)',...
                                cfnParams.numFolds,cfnParams.numRepeats);
            accuracyBit = sprintf('%.3f +/- %.3f%s',mean(meanAcc),mean(stdAcc),...
                                        cfnParams.whatLossUnits);
        else
            whatStat = sprintf('Mean +/- std (mean across %u folds; std of this mean across %u repeats)',...
                                cfnParams.numFolds,cfnParams.numRepeats);
            accuracyBit = sprintf('%.3f +/- %.3f%s',mean(meanAcc),std(meanAcc),cfnParams.whatLossUnits);
        end
    end
    fprintf(1,'\n%s %s (%u-class) using %s classification with %s:\n%s\n',...
            whatStat,cfnParams.whatLoss,cfnParams.numClasses,...
            cfnParams.whatClassifier,featureText,accuracyBit);
    if cfnParams.computePerFold
        fprintf(1,'[Training folds: mean = %.3f%s, std = %.3f%s across %u folds].\n',...
                    mean(inSampleAccStd(:,1)),cfnParams.whatLossUnits,...
                    mean(inSampleAccStd(:,2)),cfnParams.whatLossUnits,cfnParams.numFolds);
    end
    fprintf(1,'\n');
end
fprintf(1,'Looks good? Don''t forget to compare results using other (simpler) methods...\n');

%-------------------------------------------------------------------------------
%% Compute nulls for permutation testing
%-------------------------------------------------------------------------------
if numNulls > 0
    if simpleNull && (cfnParams.numFolds == 0)
        warning('The simple null is invalid for in-sample classification: using a full model-based null')
        simpleNull = false;
    end

    if simpleNull
        % Null approximates null classifer output as random label shuffles.
        % Will be a bad approximation when the classifier has e.g., class-balance bias
        % Will be a very poor approximation when CV is off (overfitting)
        fprintf(1,'Assessing %s across %u (simple permutation-based) nulls...',cfnParams.whatLoss,numNulls);
        nullStats = GiveMeSimpleNullStats(TimeSeries.Group,numNulls,cfnParams);
    else
        % Null is output from classifier from shuffled-label input.
        % More computationally intensive, but accurate.
        fprintf(1,'Computing %s across %u (full model-based) null permutations...',cfnParams.whatLoss,numNulls);

        % Compute for shuffled group label inputs:
        nullStats = zeros(numNulls,1);
        if doWithinFolds
            % The partitions of the real labels are used for every permutation, with
            % each training fold's transformation fitted once (if it does not depend
            % on the labels):
            if isfield(normSpec,'classVarFilter') && ~isempty(normSpec.classVarFilter) && normSpec.classVarFilter
                nullTransforms = {};
            else
                nullTransforms = cvOut.transforms;
            end
            nullPartitions = cvOut.partitions;
            nullLabels = cell(numNulls,1);
            for i = 1:numNulls
                nullLabels{i} = ShuffleLabels(TimeSeries.Group,nullPartitions);
            end
            if doParallel
                fprintf(1,'Using parfor for speed...\n');
                parfor i = 1:numNulls
                    nullStats(i) = BF_CrossValidate(TS_DataMat,nullLabels{i},cfnParams,normSpec,...
                                    'Operations',Operations,'partitions',nullPartitions,...
                                    'transforms',nullTransforms);
                end
            else
                for i = 1:numNulls
                    nullStats(i) = BF_CrossValidate(TS_DataMat,nullLabels{i},cfnParams,normSpec,...
                                    'Operations',Operations,'partitions',nullPartitions,...
                                    'transforms',nullTransforms);
                    if i==1
                        fprintf(1,'%u',i);
                    else
                        fprintf(1,',%u',i);
                    end
                end
            end
        elseif doParallel
            fprintf(1,'Using parfor for speed...\n');
            assert(cfnParams.computePerFold == false);
            tsGroup = TimeSeries.Group; % set up for parfor
            parfor i = 1:numNulls
                shuffledLabels = tsGroup(randperm(height(TimeSeries))); %#ok<PFBNS>
                nullStats(i) = ComputeCVAccuracies(TS_DataMat,shuffledLabels,cfnParams,true);
            end
        else
            for i = 1:numNulls
                shuffledLabels = TimeSeries.Group(randperm(height(TimeSeries)));
                nullStats(i) = ComputeCVAccuracies(TS_DataMat,shuffledLabels,cfnParams,true);
                if i==1
                    fprintf(1,'%u',i);
                else
                    fprintf(1,',%u',i);
                end
            end
        end
    end

    % Display to commandline
    fprintf(1,['\n\nMean %s (%u-class) using %u-fold %s classification with %s' ...
                     ' across %u nulls:\n%.3f +/- %.3f%%\n\n'],...
                    cfnParams.whatLoss,...
                    cfnParams.numClasses,...
                    cfnParams.numFolds,...
                    cfnParams.classifierText,...
                    featureText,...
                    numNulls,...
                    mean(nullStats),...
                    std(nullStats));

    % Plot the null distribution:
    if doPlot
        f = figure('color','w');
        hold('on')
        f.Position(3:4) = [551,277];
        h = histogram(nullStats);
        h.FaceColor = ones(1,3)*0.5;
        xline(mean(meanAcc),'r','LineWidth',2)
        if perFold
            xline(mean(inSampleAccStd(:,1)),'--r');
        end
        xlabel(sprintf('%s (%s)',cfnParams.whatLoss,cfnParams.whatLossUnits))
        ylabel('Number of nulls')
        if simpleNull
            legendEntries = {sprintf('Simple shuffled *output* label null (%u)',numNulls),...
                                        'Real labels (out-of-sample)'};
        else
            legendEntries = {sprintf('Model-based (shuffled input label) null (%u)',numNulls),...
                                    'Real labels (out-of-sample)'};
        end
        if perFold
            legendEntries{3} = 'Real labels (in-sample)';
        end
        legend(legendEntries)
    end

    % Estimate a p-value:
    % (not exactly comparing like with like but if number of nulls is high
    % enough, my intuition is that this approximation will converge to 'true'
    % value; i.e., if you did the same numRepeats for each null sample as in
    % the real sample)
    % Label-randomization null:
    % (note that for a discrete quantity like classification accuracy with a
    % small sample size, equality cases become more likely):
    pValPermTest = mean(mean(meanAcc) <= nullStats);
    fprintf(1,'Estimated p-value (permutation test) = %.2g\n',pValPermTest);

    pValZ = 1 - normcdf(mean(meanAcc),mean(nullStats),std(nullStats));
    fprintf(1,'Estimated p-value (Gaussian fit) = %.2g\n',pValZ);
else
    nullStats = [];
end

%-------------------------------------------------------------------------------
%% Plot a confusion matrix
%-------------------------------------------------------------------------------
% Convert real and predicted class labels to matrix form (numClasses x N),
% required as input to plotconfusion:
realLabels = TimeSeries.Group;

% Predict from the first CV-partition:
if doWithinFolds
    predictLabels = cvOut.yPredict(:,1);
elseif cfnParams.numFolds > 0
    % Just show results from the first CV repeat:
    predictLabels = kfoldPredict(CVMdl{1});
else
    predictLabels = predict(CVMdl{1},TS_DataMat);
end

if doPlot
    if exist('confusionchart','file') == 0 && exist('plotconfusion','file') == 0
        warning('No available confusion matrix plotting functions')
    else
        if cfnParams.numFolds > 0
            titleText = sprintf('Confusion matrix from %u-fold %s cross-validation (REPEAT 1/%u SHOWN)',...
                    cfnParams.numFolds,cfnParams.classifierText,cfnParams.numRepeats);
        else
            titleText = sprintf('Confusion matrix from in-sample %s classification',...
                                cfnParams.classifierText);
        end

        % Prefer the confusionchart to plotconfusion?
        if exist('confusionchart','file') ~= 0
            % Requires Matlab 2020 (Stats/ML Toolbox):
            figure('color','w');
            confusionchart(realLabels,predictLabels);
            title(titleText)
            % Wish I could set the interpreter to 'none', but the confusionchart
            % title text doesn't seem to have this option... :-/
        else
            % Requires the Deep Learning Toolbox:
            figure('color','w');
            plotconfusion(realLabels,predictLabels);
            % Fix axis labels:
            ax = gca;
            ax.XTickLabel(1:cfnParams.numClasses) = classLabels;
            ax.YTickLabel(1:cfnParams.numClasses) = classLabels;
            ax.TickLabelInterpreter = 'none';
            title(titleText,'interpreter','none')
        end

    end
end

%-------------------------------------------------------------------------------
%% Save trained classifier
%-------------------------------------------------------------------------------
% E.g., for loading in later for classifying a new dataset (cf. TS_Predict)
jointClassifier = []; % no classifier is trained unless a filename is specified
if doSaveClassifier
    % Train the classifier on all the data (no cross-validation):
    cfnParamsAll = cfnParams;
    cfnParamsAll.numFolds = 0;
    if doWithinFolds
        % Fit the transformation on all the data, and train on the features it keeps:
        [transformAll,XAll] = BF_FitTransform(TS_DataMat,normSpec,Operations,TimeSeries.Group);
        [inSampleAcc,Mdl] = GiveMeCfn(XAll,TimeSeries.Group,XAll,TimeSeries.Group,cfnParamsAll);
        [~,ind] = ismember(transformAll.featureParams.ID,Operations.ID);
        Operations = Operations(ind,:);
    else
        [inSampleAcc,Mdl] = GiveMeCfn(TS_DataMat,TimeSeries.Group,TS_DataMat,...
                                            TimeSeries.Group,cfnParamsAll);
    end

    jointClassifier = struct();
    jointClassifier.Mdl = Mdl;
    jointClassifier.Operation.ID = Operations.ID;
    jointClassifier.Operation.Name = Operations.Name;
    jointClassifier.classLabels = classLabels;
    jointClassifier.Accuracy = inSampleAcc;
    if cfnParams.numFolds > 0
        jointClassifier.CVAccuracy = mean(meanAcc);
    else
        jointClassifier.CVAccuracy = NaN;
    end
    if perFold
        jointClassifier.CVAccuracyFolds = reshape(allFoldData(:,2,:),cfnParams.numFolds,cfnParams.numRepeats);
    else
        jointClassifier.CVAccuracyFolds = [];
    end
    jointClassifier.whatLoss = cfnParams.whatLoss;
    jointClassifier.whatLossUnits = cfnParams.whatLossUnits;
    jointClassifier.whatClassifier = cfnParams.whatClassifier;
    jointClassifier.numFolds = cfnParams.numFolds;
    jointClassifier.numRepeats = cfnParams.numRepeats;
    if doWithinFolds
        spec = transformAll.spec;
        jointClassifier.normalizationInfo = struct('normFunction',spec.normFunction,...
                    'filterOptions',spec.filterOptions,'classVarFilter',spec.classVarFilter,...
                    'featureParams',transformAll.featureParams,'trainIDs',TimeSeries.ID,...
                    'transform',transformAll,'normalizeWithinFolds',true);
    else
        jointClassifier.normalizationInfo = ModelNormalizationInfo(TS_GetFromData(whatData,'normalizationInfo'),Operations);
    end
    classes = classLabels;

    save(cfnParams.classifierFilename,'jointClassifier','classes','cfnParams','-v7.3');
    fprintf('Saved trained classifier to ''%s''.\n',cfnParams.classifierFilename);
end

% Save my sensitive eyes an overload:
if nargout==0
    clear('meanAcc','nullStats','jointClassifier')
end

%-------------------------------------------------------------------------------
function normInfo = ModelNormalizationInfo(normInfo,Operations)
    % Normalization information for the features of the model (in its column order)
    if isempty(normInfo)
        % The data were not normalized (e.g., loaded from HCTSA.mat)
        normInfo = [];
        return
    end
    if ~isfield(normInfo,'featureParams')
        if ~ismember(normInfo.normFunction,{'none','nothing'})
            warning(['The data were normalized by a version of TS_Normalize that does not store' ...
                    ' the normalization parameters, so TS_Predict cannot use this classifier.' ...
                    ' Re-run TS_Normalize and TS_Classify to save a usable classifier.'])
        end
        return
    end
    [hasParams,ind] = ismember(Operations.ID,normInfo.featureParams.ID);
    if ~all(hasParams) || ~isequal(normInfo.featureParams.Name(ind),Operations.Name)
        warning(['The normalization parameters stored with the data do not match its features,' ...
                ' so TS_Predict cannot use this classifier.'])
        normInfo = rmfield(normInfo,'featureParams');
    else
        normInfo.featureParams = normInfo.featureParams(ind,:);
    end
    if isfield(normInfo,'transform') && isfield(normInfo,'featureParams')
        % The transformation, restricted in the same way:
        [~,ind] = ismember(Operations.ID,normInfo.transform.featureParams.ID);
        normInfo.transform.featureParams = normInfo.transform.featureParams(ind,:);
        normInfo.transform.trainMedian = normInfo.transform.trainMedian(ind);
    end
end

%-------------------------------------------------------------------------------
function [X,TimeSeries,Operations,normSpec] = LoadRawData(whatData,whatDataFile,...
                                    TimeSeries,rawData,normalizeWithinFolds)
    % The unnormalized data for the time series of whatData, and the normalization
    % specification to fit within folds
    normInfo = TS_GetFromData(whatData,'normalizationInfo');

    % The normalization specification:
    if isstruct(normalizeWithinFolds)
        normSpec = normalizeWithinFolds;
    elseif ~isempty(normInfo)
        % As recorded by TS_Normalize:
        normSpec = struct('normFunction',normInfo.normFunction,'filterOptions',normInfo.filterOptions);
        if isfield(normInfo,'classVarFilter')
            normSpec.classVarFilter = normInfo.classVarFilter;
        end
    else
        % TS_Normalize defaults:
        normSpec = struct('normFunction','mixedSigmoid','filterOptions',[0.70,1]);
    end
    if ~isfield(normSpec,'normFunction') || isempty(normSpec.normFunction)
        normSpec.normFunction = 'mixedSigmoid';
    end
    if ~isfield(normSpec,'filterOptions') || isempty(normSpec.filterOptions)
        normSpec.filterOptions = [0.70,1];
    end

    % The unnormalized data:
    if ~isempty(rawData)
        rawSource = rawData;
    elseif isempty(normInfo)
        rawSource = whatData; % the data are not normalized
    elseif ischar(whatDataFile) && endsWith(whatDataFile,'_N.mat')
        rawSource = [whatDataFile(1:end-6),'.mat'];
    else
        error(['Cannot find the unnormalized version of the data to normalize within folds:' ...
                    ' set it using the ''rawData'' option']);
    end
    [XRaw,TimeSeriesRaw,Operations,rawFile] = TS_LoadData(rawSource);
    fprintf(1,'Normalizing within folds, starting from the unnormalized data in %s\n',rawFile);

    % The time series of whatData, in the same order:
    [isFound,ind] = ismember(TimeSeries.ID,TimeSeriesRaw.ID);
    if ~all(isFound) || ~isequal(TimeSeries.Name,TimeSeriesRaw.Name(ind))
        error('The time series of %s are not all in %s (matched by ID and Name)',whatDataFile,rawFile);
    end
    X = XRaw(ind,:);
    X(~isfinite(X)) = NaN;

    % Data that were not normalized: remove time series with too few good values
    % (as TS_Normalize, using each time series' own values):
    if isempty(normInfo)
        keepRows = (1 - mean(isnan(X),2) >= normSpec.filterOptions(1));
        if ~any(keepRows)
            error('No time series had more than %4.2f%% good values',normSpec.filterOptions(1)*100);
        elseif ~all(keepRows)
            fprintf(1,'Removing %u time series with fewer than %4.2f%% good values: from %u to %u.\n',...
                        sum(~keepRows),normSpec.filterOptions(1)*100,length(keepRows),sum(keepRows));
            X = X(keepRows,:);
            TimeSeries = TimeSeries(keepRows,:);
        end
    end
end

%-------------------------------------------------------------------------------
function yShuffled = ShuffleLabels(y,partitions)
    % A random permutation of the labels with every class in every training and
    % test fold of the partitions (as for the stratified partitions of the real labels)
    numClasses = length(unique(y));
    for attempt = 1:10000
        yShuffled = y(randperm(length(y)));
        allPresent = true;
        for r = 1:length(partitions)
            for k = 1:partitions{r}.NumTestSets
                if length(unique(yShuffled(training(partitions{r},k)))) < numClasses || ...
                        length(unique(yShuffled(test(partitions{r},k)))) < numClasses
                    allPresent = false;
                end
            end
        end
        if allPresent
            return
        end
    end
    error('Could not permute the labels with every class in every training and test fold');
end

%-------------------------------------------------------------------------------
function n = TestFoldSizes(cvOut,numObs)
    % Number of observations in each test fold (numFolds x numRepeats)
    if isempty(cvOut.partitions)
        n = numObs;
    else
        n = cell2mat(cellfun(@(c) c.TestSize(:),cvOut.partitions(:)','UniformOutput',false));
    end
end

end
