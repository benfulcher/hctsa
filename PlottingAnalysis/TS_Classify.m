function [foldLosses,nullStat,jointClassifier] = TS_Classify(whatData,cfnParams,varargin)
% TS_Classify   Classify groups in the data using all features
%
% TS_Classify uses a classifier to learn group labels assigned to time series
% in the dataset.
%
% It can also provide insights into the dimensionality required for good performance
% by classifying with the leading principal components of the feature space.
%
%---USAGE:
% TS_Classify();
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData)
% cfnParams, parameters for the classification algorithm (cf GiveMeDefaultClassificationParams)
%
% (OPTIONAL):
% 'numNulls' (numeric), number of nulls to compute (0 for no null comparison)
% 'seedReset', input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)
% 'doPlot', whether to plot outputs (including confusion matrix)
%
%---OUTPUTS:
% Text output on classification rate using all features, and if doPCs = true, also
% shows how this varies as a function of reduced PCs (text and as an output plot)
% 'foldLosses', the performance metric across repeats of cross-validation
% 'nullStat', the performance metric across randomizations of the data labels

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
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
% <set default classificaiton parameters after inspecting TimeSeries labeling below>

% Use an inputParser to parse additional options as parameters:
inputP = inputParser;
% numNulls:
default_numNulls = 20;
check_numNulls = @(x) isnumeric(x);
addParameter(inputP,'numNulls',default_numNulls,check_numNulls);
% seedReset:
default_seedReset = 'default';
check_seedReset = @(x) ischar(x);
addParameter(inputP,'seedReset',default_seedReset,check_seedReset);
% doPlot:
default_doPlot = true;
check_doPlot = @(x) islogical(x);
addParameter(inputP,'doPlot',default_doPlot,check_doPlot);

% Parse input arguments:
parse(inputP,varargin{:});
numNulls = inputP.Results.numNulls;
seedReset = inputP.Results.seedReset;
doPlot = inputP.Results.doPlot;
clear('inputP');

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series in %s. Use TS_LabelGroups.',whatDataFile);
end

if nargin < 2
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end
groupNames = TS_GetFromData(whatData,'groupNames');
numFeatures = height(Operations);

%-------------------------------------------------------------------------------
% Give basic info about the represented classes:
dataStruct = makeDataStruct();
TellMeAboutLabeling(dataStruct);

%-------------------------------------------------------------------------------
% Fit the model
%-------------------------------------------------------------------------------
% Reset the random seed:
BF_ResetSeed(seedReset); % reset the random seed for CV-reproducibility

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance
TellMeAboutClassification(cfnParams);

CVMdl = cell(cfnParams.numRepeats,1);
foldLosses = zeros(cfnParams.numRepeats,1);
cfnParams.computePerFold = false;
for i = 1:cfnParams.numRepeats
    [foldLosses(i),CVMdl{i}] = GiveMeFoldLosses(TS_DataMat,TimeSeries.Group);
end

% Display results to commandline:
if cfnParams.numRepeats==1
    foldBit = sprintf('%u-fold',cfnParams.numFolds);
    accuracyBit = sprintf('%.3f%%',mean(foldLosses));
else
    foldBit = sprintf('%u-fold (%u repeats)',cfnParams.numFolds,cfnParams.numRepeats);
    accuracyBit = sprintf('%.3f +/- %.3f%%',mean(foldLosses),std(foldLosses));
end
fprintf(1,['\nMean (across folds) %s (%u-class) using %s %s classification with %u' ...
                 ' features:\n%s\n\n'],...
        cfnParams.whatLoss,cfnParams.numClasses,foldBit,...
        cfnParams.whatClassifier,numFeatures,accuracyBit);

% Plot as a histogram:
if doPlot
    f = figure('color','w');
    histogram(foldLosses*100)
    xlim([0,100]);
end

%-------------------------------------------------------------------------------
% Check nulls:
%-------------------------------------------------------------------------------
if numNulls > 0
    fprintf(1,'Computing %s across %u null permutations...',cfnParams.whatLoss,numNulls);
    nullStat = zeros(numNulls,1);
    parfor i = 1:numNulls
        % Compute for shuffled data labels:
        shuffledLabels = TimeSeries.Group(randperm(height(TimeSeries)));
        nullStat(i) = GiveMeCfn(TS_DataMat,shuffledLabels,[],[],cfnParams);
        if i==1
            fprintf(1,'%u',i);
        else
            fprintf(1,',%u',i);
        end
    end
    fprintf(1,['\n\nMean %s (%u-class) using %u-fold %s classification with %u' ...
                     ' features across %u nulls:\n%.3f +/- %.3f%%\n\n'],...
                    cfnParams.whatLoss,...
                    cfnParams.numClasses,...
                    cfnParams.numFolds,...
                    cfnParams.whatClassifier,...
                    numFeatures,...
                    numNulls,...
                    mean(nullStat),...
                    std(nullStat));

    % Plot the null distribution:
    if doPlot
        f = figure('color','w'); hold('on')
        ax = gca;
        histogram(nullStat,'normalization','pdf');
        plot(ones(2,1)*mean(foldLosses),ax.YLim,'r','LineWidth',2)
        xlabel(cfnParams.whatLoss)
        ylabel('Probability density')
        legend('shuffled labels','real labels')
    end

    % Estimate a p-value:
    % (not exactly comparing like with like but if number of nulls is high
    % enough, my intuition is that this approximation will converge to 'true'
    % value; i.e., if you did the same numRepeats for each null sample as in
    % the real sample)
    % Label-randomization null:
    % (note that for a discrete quantity like classification accuracy with a
    % small sample size, equality cases become more likely):
    pValPermTest = mean(mean(foldLosses) <= nullStat);
    fprintf(1,'Estimated p-value (permutation test) = %.2g\n',pValPermTest);

    pValZ = 1 - normcdf(mean(foldLosses),mean(nullStat),std(nullStat));
    fprintf(1,'Estimated p-value (Gaussian fit) = %.2g\n',pValZ);
else
    nullStat = [];
end

%-------------------------------------------------------------------------------
% Plot a confusion matrix
%-------------------------------------------------------------------------------
% Convert real and predicted class labels to matrix form (numClasses x N),
% required as input to plotconfusion:
realLabels = BF_ToBinaryClass(TimeSeries.Group,numClasses,false);
% Predict from the first CV-partition:
predictLabels = BF_ToBinaryClass(kfoldPredict(CVMdl{1}),cfnParams.numClasses,false);
if doPlot
    try
        if exist('confusionchart','file') ~= 0
            % Requires Matlab 2020 (Stats/ML Toolbox):
            f = figure('color','w');
            confusionchart(realLabels,predictLabels);
        else
            % Requires the Deep Learning Toolbox:
            plotconfusion(realLabels,predictLabels);
        end

        % Fix axis labels:
        ax = gca;
        ax.XTickLabel(1:numClasses) = groupNames;
        ax.YTickLabel(1:numClasses) = groupNames;
        ax.TickLabelInterpreter = 'none';

        title(sprintf('Confusion matrix from 10-fold (%u repeats) cross-validation',numRepeats));
    catch
        warning('No available confusion matrix plotting functions')
    end
end

%-------------------------------------------------------------------------------
% Save classifier:
%-------------------------------------------------------------------------------
if ~isempty(cfnParams.classifierFilename)
    [Acc,Mdl] = GiveMeCfn(TS_DataMat,TimeSeries.Group,TS_DataMat,...
                                        TimeSeries.Group,cfnParams);

    Operations = TS_GetFromData(whatData,'Operations');
    jointClassifier.Operation.ID = Operations.ID;
    jointClassifier.Operation.Name = Operations.Name;
    jointClassifier.CVAccuracy = mean(foldLosses); % Cross-validated accuracy
    jointClassifier.Accuracy = Acc; % For some reason the Acc retrieved as above is amazing?
    jointClassifier.Mdl = Mdl;
    jointClassifier.whatLoss = whatLoss;
    jointClassifier.normalizationInfo = TS_GetFromData(whatData,'normalizationInfo');
    classes = groupNames;

    if exist(classifierFilename,'file')
        out = input(sprintf(['File %s already exists -- continuing will overwrite the file.'...
                  '\n[Press ''y'' to continue] '],classifierFilename), 's');
        if ~strcmp(out,'y')
            return;
        end
    end
    save(classifierFilename,'jointClassifier','classes','cfnParams','-v7.3');
    fprintf('Saved trained classifier to ''%s''.\n',classifierFilename);
end

%-------------------------------------------------------------------------------
function [foldLosses,CVMdl] = GiveMeFoldLosses(dataMatrix,dataLabels)
    % Returns the output (e.g., loss) for the custom fn_loss function across all folds
    [foldLosses,CVMdl] = GiveMeCfn(dataMatrix,dataLabels,[],[],cfnParams);
end
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
