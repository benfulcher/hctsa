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
% numNulls (numeric), number of nulls to compute (0 for no null comparison)
%
% (OPTIONAL):
% doParallel [false], whether to speed up the null computation using a parfor loop
% doPlot [true], whether to plot outputs (including confusion matrix)
% seedReset, input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)

%
%---OUTPUTS:
% Text output on classification rate using all features, and if doPCs = true, also
% shows how this varies as a function of reduced PCs (text and as an output plot)
% foldLosses, the performance metric across repeats of cross-validation
% nullStats, the performance metric across randomizations of the data labels
% jointClassifier, details of the saved all-features classifier

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
% <set default classification parameters after inspecting TimeSeries labeling below>
if nargin < 3
    numNulls = 0;
end

% Use an inputParser to parse additional options as parameters:

inputP = inputParser;
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

% Parse input arguments:
parse(inputP,varargin{:});
doParallel = inputP.Results.doParallel;
seedReset = inputP.Results.seedReset;
doPlot = inputP.Results.doPlot;
clear('inputP');

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
[groupLabels,classLabels,groupLabelsInteger,numGroups] = TS_ExtractGroupLabels(TimeSeries);

% Give basic info about the represented classes:
TellMeAboutLabeling(TimeSeries);

% Settings for the classification model:
if nargin < 2 || (isstruct(cfnParams) && isempty(fields(cfnParams))) || isempty(cfnParams)
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end
TellMeAboutClassification(cfnParams);

% Filter down a reduced feature set if required:
[TS_DataMat,Operations] = FilterFeatures(TS_DataMat,Operations,cfnParams);
numFeatures = height(Operations);

%-------------------------------------------------------------------------------
% Fit the model
%-------------------------------------------------------------------------------
% Reset the random seed for CV-reproducibility
BF_ResetSeed(seedReset);

% Fit the classification model to the dataset and evaluate performance:
CVMdl = cell(cfnParams.numRepeats,1);
foldLosses = zeros(cfnParams.numRepeats,1);
cfnParams.computePerFold = false;
for i = 1:cfnParams.numRepeats
    [meanAcc(i),CVMdl{i}] = GiveMeFoldLosses(TS_DataMat,TimeSeries.Group);
end

% Display results to commandline:
if cfnParams.numFolds == 0
    assert(cfnParams.numRepeats==1)
    fprintf(1,'In-sample %s using %s = %.3f%s\n',cfnParams.whatLoss,cfnParams.whatClassifier,...
                                            meanAcc,cfnParams.whatLossUnits);
else
    if cfnParams.numRepeats==1
        whatStat = sprintf('Mean (across %u folds)',cfnParams.numFolds);
        accuracyBit = sprintf('%.3f%s',meanAcc,cfnParams.whatLossUnits);
    else
        whatStat = sprintf('Mean (across %u folds, then across %u repeats)',cfnParams.numFolds,cfnParams.numRepeats);
        accuracyBit = sprintf('%.3f +/- %.3f%s',mean(meanAcc),std(meanAcc),cfnParams.whatLossUnits);
    end
    fprintf(1,['\n%s %s (%u-class) using %s classification with %u' ...
                     ' features:\n%s\n\n'],...
            whatStat,cfnParams.whatLoss,cfnParams.numClasses,...
            cfnParams.whatClassifier,numFeatures,accuracyBit);
end


%-------------------------------------------------------------------------------
% Check nulls:
%-------------------------------------------------------------------------------
if numNulls > 0
    fprintf(1,'Computing %s across %u null permutations...',cfnParams.whatLoss,numNulls);
    nullStats = zeros(numNulls,1);
    % Compute for shuffled data labels:
    if doParallel
        fprintf(1,'Using parfor for speed...\n');
        parfor i = 1:numNulls
            shuffledLabels = TimeSeries.Group(randperm(height(TimeSeries)));
            nullStats(i) = GiveMeCfn(TS_DataMat,shuffledLabels,TS_DataMat,shuffledLabels,cfnParams);
        end
    else
        for i = 1:numNulls
            shuffledLabels = TimeSeries.Group(randperm(height(TimeSeries)));
            nullStats(i) = GiveMeCfn(TS_DataMat,shuffledLabels,TS_DataMat,shuffledLabels,cfnParams);
            if i==1
                fprintf(1,'%u',i);
            else
                fprintf(1,',%u',i);
            end
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
                    mean(nullStats),...
                    std(nullStats));

    % Plot the null distribution:
    if doPlot
        f = figure('color','w');
        hold('on')
        f.Position(3:4) = [551,277];
        ax = gca;
        h = histogram(nullStats);
        h.FaceColor = ones(1,3)*0.5;
        plot(ones(2,1)*mean(foldLosses),ax.YLim,'r','LineWidth',2)
        xlabel(sprintf('%s (%s)',cfnParams.whatLoss,cfnParams.whatLossUnits))
        ylabel('Number of nulls')
        legend(sprintf('Shuffled labels (%u)',numNulls),'Real labels')
    end

    % Estimate a p-value:
    % (not exactly comparing like with like but if number of nulls is high
    % enough, my intuition is that this approximation will converge to 'true'
    % value; i.e., if you did the same numRepeats for each null sample as in
    % the real sample)
    % Label-randomization null:
    % (note that for a discrete quantity like classification accuracy with a
    % small sample size, equality cases become more likely):
    pValPermTest = mean(mean(foldLosses) <= nullStats);
    fprintf(1,'Estimated p-value (permutation test) = %.2g\n',pValPermTest);

    pValZ = 1 - normcdf(mean(foldLosses),mean(nullStats),std(nullStats));
    fprintf(1,'Estimated p-value (Gaussian fit) = %.2g\n',pValZ);
else
    nullStats = [];
end

%-------------------------------------------------------------------------------
% Plot a confusion matrix
%-------------------------------------------------------------------------------
% Convert real and predicted class labels to matrix form (numClasses x N),
% required as input to plotconfusion:
realLabels = TimeSeries.Group;
% Predict from the first CV-partition:
if cfnParams.numFolds > 0
    % Just show results from the first CV repeat:
    predictLabels = kfoldPredict(CVMdl{1});
else
    predictLabels = predict(CVMdl{1},TS_DataMat);
end

if doPlot
    if exist('confusionchart','file') == 0 && exist('plotconfusion','file') ==0
        warning('No available confusion matrix plotting functions')
    else
        if cfnParams.numFolds > 0
            titleText = sprintf('Confusion matrix from %u-fold %s cross-validation (repeat 1)',...
                    cfnParams.numFolds,cfnParams.whatClassifier);
        else
            titleText = sprintf('Confusion matrix from in-sample %s classification',...
                                cfnParams.whatClassifier);
        end

        % Prefer the confusionchart to plotconfusion?
        if exist('confusionchart','file') ~= 0
            % Requires Matlab 2020 (Stats/ML Toolbox):
            f_confusion = figure('color','w');
            confusionchart(realLabels,predictLabels);
            title(titleText)
            % Wish I could set the interpreter to 'none', but the confusionchart
            % title text doesn't seem to have this option... :-/
        else
            % Requires the Deep Learning Toolbox:
            f_confusion = figure('color','w');
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
% Save trained classifier:
%-------------------------------------------------------------------------------
if ~isempty(cfnParams.classifierFilename)
    [Acc,Mdl] = GiveMeCfn(TS_DataMat,TimeSeries.Group,TS_DataMat,...
                                            TimeSeries.Group,cfnParams);
    jointClassifier.Operation.ID = Operations.ID;
    jointClassifier.Operation.Name = Operations.Name;
    jointClassifier.CVAccuracy = mean(foldLosses);
    jointClassifier.Accuracy = Acc;
    jointClassifier.Mdl = Mdl;
    jointClassifier.whatLoss = whatLoss;
    jointClassifier.normalizationInfo = TS_GetFromData(whatData,'normalizationInfo');
    classes = classLabels;

    if exist(cfnParams.classifierFilename,'file')
        out = input(sprintf(['File %s already exists -- continuing will overwrite the file.'...
                  '\n[Press ''y'' to continue] '],cfnParams.classifierFilename), 's');
        if ~strcmp(out,'y')
            return;
        end
    end
    save(cfnParams.classifierFilename,'jointClassifier','classes','cfnParams','-v7.3');
    fprintf('Saved trained classifier to ''%s''.\n',cfnParams.classifierFilename);
end

% Save my sensitive eyes an overload:
if nargout==0
    clear('foldLosses','nullStats','jointClassifier')
end

%-------------------------------------------------------------------------------
function [foldLosses,CVMdl] = GiveMeFoldLosses(dataMatrix,dataLabels)
    % Returns the output (e.g., loss) for the custom fn_loss function across all folds
    [foldLosses,CVMdl] = GiveMeCfn(dataMatrix,dataLabels,dataMatrix,dataLabels,cfnParams);
end
%-------------------------------------------------------------------------------

end
