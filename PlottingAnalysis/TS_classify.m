function TS_classify(whatData,whatClassifier,doPCs,seedReset)
% TS_classify   Classify groups in the data using all features (and PCs)
%
% This function uses a classifier to learn group labels assigned to time series
% in the dataset.
%
%---USAGE:
% TS_classify();
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData)
% whatClassifier, the classification method to use (e.g., 'svm_linear', 'knn', 'linear')
% doPCs, (binary) whether to investigate the behavior of different numbers of PCs of the
%        data matrix (default: 1)
% seedReset, input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)
%
%---OUTPUTS:
% Text output on classification rate using all features, and if doPCs = 1, also
% shows how this varies as a function of reduced PCs (text and as an output plot)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
    whatClassifier = 'svm_linear';
    % 'svm', 'discriminant', 'knn'
end
if nargin < 3
    doPCs = 1;
end
if nargin < 4
    seedReset = 'default';
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,dataFile] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~isfield(TimeSeries,'Group')
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
groupNames = TS_GetFromData(whatData,'groupNames');
timeSeriesGroup = [TimeSeries.Group]'; % Use group form (column vector)
numClasses = max(timeSeriesGroup); % assuming group in form of integer class labels starting at 1
numFeatures = length(Operations);
numTimeSeries = length(TimeSeries);

%-------------------------------------------------------------------------------
% Fit the model
%-------------------------------------------------------------------------------
numFolds = howManyFolds(timeSeriesGroup,numClasses);

BF_ResetSeed(seedReset); % reset the random seed for CV-reproducibility

%-------------------------------------------------------------------------------
% Fit the classification model to the dataset (for each cross-validation fold)
% and evaluate performance
fprintf(1,['Training and evaluating a %u-class %s classifier in a %u-feature' ...
                        ' space using %u-fold cross validation...\n'],...
                        numClasses,whatClassifier,numFeatures,numFolds);
[foldLosses,CVMdl,outputStat] = GiveMeFoldLosses(TS_DataMat,timeSeriesGroup);

fprintf(1,['\n%s (%u-class) using %u-fold %s classification with %u' ...
                 ' features:\n%.3f +/- %.3f%%\n\n'],...
                    outputStat,...
                    numClasses,...
                    numFolds,...
                    whatClassifier,...
                    numFeatures,...
                    mean(foldLosses),...
                    std(foldLosses));

% f = figure('color','w');
% histogram(foldLosses*100)
% xlim([0,100]);

%-------------------------------------------------------------------------------
% Check nulls:
%-------------------------------------------------------------------------------
doNull = 0;
numNulls = 20;
if doNull
    meanNull = zeros(numNulls,1);
    for i = 1:numNulls
        % Compute for shuffled data labels:
        shuffledLabels = timeSeriesGroup(randperm(length(timeSeriesGroup)));
        foldLosses_null = GiveMeFoldLosses(TS_DataMat,shuffledLabels);
        meanNull(i) = mean(foldLosses_null);
    end
    fprintf(1,['\nMean %s (%u-class) using %u-fold %s classification with %u' ...
                     ' features across %u nulls:\n%.3f +/- %.3f%%\n\n'],...
                    outputStat,...
                    numClasses,...
                    numFolds,...
                    whatClassifier,...
                    numFeatures,...
                    numNulls,...
                    mean(meanNull),...
                    std(meanNull));
end

%-------------------------------------------------------------------------------
% Plot confusion matrix
%-------------------------------------------------------------------------------
% CONVERT BOTH TO numClassesxN form before plotting confusion matrix
realLabels = BF_ToBinaryClass(timeSeriesGroup');
predictLabels = BF_ToBinaryClass(kfoldPredict(CVMdl));
plotconfusion(realLabels,predictLabels);

% Fix axis labels:
ax = gca;
ax.XTickLabel(1:numClasses) = groupNames;
ax.YTickLabel(1:numClasses) = groupNames;
ax.TickLabelInterpreter = 'none';
% Make a nice white figure background:
f = gcf; f.Color = 'w';

%-------------------------------------------------------------------------------
% Compare performance of reduced PCs from the data matrix:
%-------------------------------------------------------------------------------
if doPCs
    if any(isnan(TS_DataMat(:)))
        warning(['Cannot compute PCs of data matrix containing NaNs...\n' ...
                    '(to compute PCs, re-run TS_Normalize to filter out all NaNs)'])
        return
    end

    % Compute top 10 PCs of the data matrix:
    numPCs = 10;
    fprintf('Computing top %u PCs...',numPCs)
    [~, pcScore, ~, ~, ~] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
    fprintf(' Done.\n')
    numPCs = min(10,size(pcScore,2)); % sometimes lower than attempted 10

    % Compute cumulative performance of PCs:
    cfnRate = zeros(numPCs,2);
    fprintf('Computing classification rates keeping top 1-%u PCs...\n',numPCs)
    for i = 1:numPCs
        PCfoldLosses = GiveMeFoldLosses(pcScore(:,1:i),timeSeriesGroup);
        cfnRate(i,1) = mean(PCfoldLosses);
        cfnRate(i,2) = std(PCfoldLosses);
        fprintf(1,'%u PCs:   %.3f +/- %.3f%%\n',i,cfnRate(i,1),cfnRate(i,2));
    end

    plotColors = BF_getcmap('spectral',3,1);

    f = figure('color','w'); hold on
    plot([1,numPCs],ones(2,1)*mean(foldLosses),'--','color',plotColors{3})
    plot(1:numPCs,cfnRate(:,1),'o-k')
    legend(sprintf('All %u features (%.1f%%)',numFeatures,mean(foldLosses)),...
                sprintf('PCs (%.1f--%.1f%%)',min(cfnRate(:,1)),max(cfnRate(:,1))))
    plot(1:numPCs,cfnRate(:,1)+cfnRate(:,2),':k')
    plot(1:numPCs,cfnRate(:,1)-cfnRate(:,2),':k')

    xlabel('Number of PCs');
    ylabel('Classification accuracy (%)')

    titleText = sprintf('Classification rate (%u-class) using %u-fold %s classification',...
                                numClasses,...
                                numFolds,...
                                whatClassifier);
    title(titleText)

end

%-------------------------------------------------------------------------------
function [foldLosses,CVMdl,whatLoss] = GiveMeFoldLosses(dataMatrix,dataLabels)
    % Returns the output (e.g., loss) for the custom fn_loss function across all folds
    [foldLosses,CVMdl,whatLoss] = GiveMeCfn(whatClassifier,dataMatrix,...
                            dataLabels,[],[],numClasses,[],[],1,numFolds);
end

end
