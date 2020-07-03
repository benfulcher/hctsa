function [foldLosses,nullStat,jointClassifier] = TS_Classify(whatData,whatClassifier,varargin)
% TS_Classify   Classify groups in the data using all features
%
% This function uses a classifier to learn group labels assigned to time series
% in the dataset.
% It can also provide insights into the dimensionality required for good performance
% by classifying with the leading principal components of the feature space
%
%---USAGE:
% TS_Classify();
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData)
% whatClassifier, the classification method to use (e.g., 'svm_linear', 'knn', 'linear')
%
% (OPTIONAL):
% 'numPCs', investigate classification using up to this many PCs of the data
%              matrix (default: 0).
% 'numNulls' (numeric), number of nulls to compute (0 for no null comparison)
% 'seedReset', input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)
% 'numFolds', number of cross-validation folds to use.
% 'numRepeats', number of times to repeat the prediction across different data partitions.
% 'classifierFilename', .mat file to save the classifier to (not saved if empty).
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
if nargin < 2
    whatClassifier = 'svm_linear'; % 'svm', 'discriminant', 'knn'
end

% Use an inputParser to control additional plotting options as parameters:
inputP = inputParser;
% numPCs:
default_numPCs = 0;
check_numPCs = @(x) isnumeric(x);
addParameter(inputP,'numPCs',default_numPCs,check_numPCs);
% numNulls:
default_numNulls = 20;
check_numNulls = @(x) isnumeric(x);
addParameter(inputP,'numNulls',default_numNulls,check_numNulls);
% seedReset:
default_seedReset = 'default';
check_seedReset = @(x) ischar(x);
addParameter(inputP,'seedReset',default_seedReset,check_seedReset);
% numFolds:
default_numFolds = [];
check_numFolds = @(x) isnumeric(x);
addParameter(inputP,'numFolds',default_numFolds,check_numFolds);
% numRepeats:
default_numRepeats = 1;
check_numRepeats = @(x) isnumeric(x);
addParameter(inputP,'numRepeats',default_numRepeats,check_numRepeats);
% doPlot:
default_doPlot = true;
check_doPlot = @(x) islogical(x);
addParameter(inputP,'doPlot',default_doPlot,check_doPlot);
% classifierFilename (save trained classifer out to file):
default_classifierFilename = '';
check_classifierFilename = @(x) ischar(x);
addParameter(inputP,'classifierFilename',default_classifierFilename,check_classifierFilename);

% Parse input arguments:
parse(inputP,varargin{:});
numPCs = inputP.Results.numPCs;
numNulls = inputP.Results.numNulls;
seedReset = inputP.Results.seedReset;
numFolds = inputP.Results.numFolds;
numRepeats = inputP.Results.numRepeats;
doPlot = inputP.Results.doPlot;
classifierFilename = inputP.Results.classifierFilename;
clear('inputP');

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
groupNames = TS_GetFromData(whatData,'groupNames');
if any(TimeSeries.Group==0)
    error('Error labeling time series groups');
end
numClasses = max(TimeSeries.Group); % assuming group in form of integer class labels starting at 1
numFeatures = height(Operations);

%-------------------------------------------------------------------------------
% Give basic info about the represented classes:
fprintf(1,'%u-class classification using %s:\n',numClasses,whatClassifier);
for i = 1:numClasses
    fprintf(1,'%u time series of class ''%s''\n',sum(TimeSeries.Group==i),groupNames{i});
end

%-------------------------------------------------------------------------------
% Fit the model
%-------------------------------------------------------------------------------
if isempty(numFolds) || numFolds==0
    % Use a heuristic to set a default number of folds given the data set size,
    % number of classes
    numFolds = HowManyFolds(TimeSeries.Group,numClasses);
end

% Reset the random seed:
BF_ResetSeed(seedReset); % reset the random seed for CV-reproducibility

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance
fprintf(1,['Training and evaluating a %u-class %s classifier in a %u-feature' ...
                ' space using %u-fold cross validation (with %u repeat(s))...\n'],...
                numClasses,whatClassifier,numFeatures,numFolds,numRepeats);
CVMdl = cell(numRepeats,1);
foldLosses = zeros(numRepeats,1);
for i = 1:numRepeats
    [foldLosses(i),CVMdl{i},outputStat] = GiveMeFoldLosses(TS_DataMat,TimeSeries.Group);
end

% Display results to commandline:
if numRepeats==1
    foldBit = sprintf('%u-fold',numFolds);
    accuracyBit = sprintf('%.3f%%',mean(foldLosses));
else
    foldBit = sprintf('%u-fold (%u repeats)',numFolds,numRepeats);
    accuracyBit = sprintf('%.3f +/- %.3f%%',mean(foldLosses),std(foldLosses));
end
fprintf(1,['\nMean (across folds) %s (%u-class) using %s %s classification with %u' ...
                 ' features:\n%s\n\n'],...
        outputStat,numClasses,foldBit,whatClassifier,numFeatures,accuracyBit);

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
    fprintf(1,'Computing %s across %u null permutations...',outputStat,numNulls);
    nullStat = zeros(numNulls,1);
    parfor i = 1:numNulls
        % Compute for shuffled data labels:
        shuffledLabels = TimeSeries.Group(randperm(height(TimeSeries)));
        nullStat(i) = GiveMeCfn(whatClassifier,TS_DataMat,shuffledLabels,[],[],numClasses,[],[],true,numFolds);
        % nullStat(i) = GiveMeFoldLosses(TS_DataMat,shuffledLabels);
        % meanNull(i) = mean(foldLosses_null);
        if i==1
            fprintf(1,'%u',i);
        else
            fprintf(1,',%u',i);
        end
    end
    fprintf(1,['\n\nMean %s (%u-class) using %u-fold %s classification with %u' ...
                     ' features across %u nulls:\n%.3f +/- %.3f%%\n\n'],...
                    outputStat,...
                    numClasses,...
                    numFolds,...
                    whatClassifier,...
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
        xlabel(outputStat)
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
predictLabels = BF_ToBinaryClass(kfoldPredict(CVMdl{1}),numClasses,false);
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
% Compare performance of reduced PCs from the data matrix:
%-------------------------------------------------------------------------------
if numPCs > 0
    if any(isnan(TS_DataMat(:)))
        warning(['Cannot compute PCs of data matrix containing NaNs...\n' ...
                    '(to compute PCs, re-run TS_Normalize to filter out all NaNs)'])
        return
    end

    % Compute top X PCs of the data matrix:
    fprintf('Computing top %u PCs...',numPCs)
    [~,pcScore,~,~,~] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
    fprintf(' Done.\n')
    numPCs = min(numPCs,size(pcScore,2)); % sometimes lower than attempted 10

    % Compute cumulative performance of PCs:
    cfnRate = zeros(numPCs,1);
    fprintf('Computing classification rates keeping top 1-%u PCs...\n',numPCs)
    for i = 1:numPCs
        cfnRate(i) = GiveMeFoldLosses(pcScore(:,1:i),TimeSeries.Group);
        fprintf(1,'%u PCs:  %.3f%%\n',i,cfnRate(i));
        % fprintf(1,'%u PCs:   %.3f +/- %.3f%%\n',i,cfnRate(i,1),cfnRate(i,2));
    end

    if doPlot
        plotColors = BF_GetColorMap('spectral',3,1);
        f = figure('color','w'); hold on
        plot([1,numPCs],ones(2,1)*mean(foldLosses),'--','color',plotColors{3})
        plot(1:numPCs,cfnRate(:,1),'o-k')
        legend(sprintf('All %u features (%.1f%%)',numFeatures,mean(foldLosses)),...
                    sprintf('PCs (%.1f-%.1f%%)',min(cfnRate),max(cfnRate)))
        % plot(1:numPCs,cfnRate(:,1)+cfnRate(:,2),':k')
        % plot(1:numPCs,cfnRate(:,1)-cfnRate(:,2),':k')

        xlabel('Number of PCs');
        ylabel('Classification accuracy (%)')

        titleText = sprintf('Classification rate (%u-class) using %u-fold %s classification',...
                                    numClasses,numFolds,whatClassifier);
        title(titleText,'interpreter','none')
    end
end

%-------------------------------------------------------------------------------
% Save classifier:
%-------------------------------------------------------------------------------
if ~isempty(classifierFilename)
    [Acc,Mdl,whatLoss] = GiveMeCfn(whatClassifier,TS_DataMat,...
                            TimeSeries.Group,TS_DataMat,TimeSeries.Group,...
                            numClasses,[],[],true,0);

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
                  '\n[Press ''y'' to continue] '], classifierFilename), 's');
        if ~strcmp(out,'y')
            return;
        end
    end
    save(classifierFilename,'jointClassifier','classes','-v7.3');
    fprintf('Saved trained classifier to ''%s''.\n',classifierFilename);
end

%-------------------------------------------------------------------------------
function [foldLosses,CVMdl,whatLoss] = GiveMeFoldLosses(dataMatrix,dataLabels)
    % Returns the output (e.g., loss) for the custom fn_loss function across all folds
    [foldLosses,CVMdl,whatLoss] = GiveMeCfn(whatClassifier,dataMatrix,...
                            dataLabels,[],[],numClasses,[],[],true,numFolds);
end

end
