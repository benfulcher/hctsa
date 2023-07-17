function TS_ClassifyLowDim(whatData,cfnParams,numPCs,doComputeAtMaxInSample)
% TS_ClassifyLowDim compare performance of reduced PCs from the data matrix
%-------------------------------------------------------------------------------
% 'whatData', HCTSA data file (or structure).
% 'cfnParams', parameters of the classification to be performed (cf.
%                   GiveMeDefaultClassificationParams)
% 'numPCs', investigate classification using up to this many PCs of the data
%              matrix (default: 5).
% 'docomputeAtMaxInsample', compute accuracy at a a number of PCs that hits 100%
%                           in-sample accuracy.

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

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
if nargin < 3
    numPCs = 10;
end
if nargin < 4
    doComputeAtMaxInSample = true;
end

%-------------------------------------------------------------------------------
%% Load data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
[groupLabels,classLabels,groupLabelsInteger,numGroups] = TS_ExtractGroupLabels(TimeSeries);
TellMeAboutLabeling(TimeSeries);

if nargin < 2 || isempty(cfnParams)
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end
TellMeAboutClassification(cfnParams);

% Filter down to reduced features if specified/required:
[TS_DataMat,Operations] = FilterFeatures(TS_DataMat,Operations,cfnParams);
numFeatures = size(TS_DataMat,2);

if doComputeAtMaxInSample
    if ~cfnParams.computePerFold
        cfnParams.computePerFold = true;
        warning(['You must compute in-sample accuracies (per fold) for',...
                        'doComputeAtMaxInSample to work: setting to true.']);
    end
end

%-------------------------------------------------------------------------------
%% Check for NaNs in data matrix
%-------------------------------------------------------------------------------
if any(isnan(TS_DataMat(:)))
    warning('Consider re-running TS_Normalize to filter out NaNs from the feature matrix)')
end

%-------------------------------------------------------------------------------
%% Compute leading X PCs of the data matrix
%-------------------------------------------------------------------------------
% Use pca to compute the first two principal components:
% (project data into space of PC scores, Y)
if ~any(isnan(TS_DataMat))
    fprintf('Computing top %u PCs...',numPCs)
    if doComputeAtMaxInSample
        [pcCoeff,pcScore,~,~,percVar] = pca(zscore(TS_DataMat));
    else
        [pcCoeff,pcScore,~,~,percVar] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
    end
else
    warning(['Data matrix contains %.2g%% NaNs. Estimating covariances on remaining data...\n' ...
                '(Could take some time...)'],100*mean(isnan(TS_DataMat(:))))
    if doComputeAtMaxInSample
        [pcCoeff,pcScore,~,~,percVar] = pca(BF_NormalizeMatrix(TS_DataMat,'zscore'),...
                            'Rows','pairwise');
    else
        [pcCoeff,pcScore,~,~,percVar] = pca(BF_NormalizeMatrix(TS_DataMat,'zscore'),...
                            'Rows','pairwise','NumComponents',numPCs);
    end
    % If this fails (covariance matrix not positive definite), can try the
    % (...,'algorithm','als') option in pca... (or toolbox for probabilistic PCA)
end
fprintf(' Done.\n')
numPCs = min(numPCs,size(pcScore,2)); % sometimes lower than the number attempted

%-------------------------------------------------------------------------------
%% Display some info about feature loading onto the reduced components
%-------------------------------------------------------------------------------
numTopLoadFeat = min(numFeatures,20); % display this many features loading onto each PC
LowDimDisplayTopLoadings(numTopLoadFeat,numPCs,pcCoeff,pcScore,TS_DataMat,Operations);

%-------------------------------------------------------------------------------
%% Compute performance with leading PCs of the feature space:
%-------------------------------------------------------------------------------
cfnRatePCs = zeros(numPCs,1);
stdAcc = zeros(numPCs,1);
if cfnParams.computePerFold
    inSampleStats = zeros(numPCs,2);
    allFoldPCs = cell(numPCs,1);
end
fprintf('Computing %s, keeping top 1-%u PCs...\n',cfnParams.whatLoss,numPCs)
cfnParams.suppressWarning = true;
for i = 1:numPCs
    topPCs = pcScore(:,1:i);
    if cfnParams.computePerFold
        [cfnRatePCs(i),stdAcc(i),inSampleStats(i,:),allFoldPCs{i}] = ...
                        ComputeCVAccuracies(topPCs,TimeSeries.Group,cfnParams,true);
        fprintf(1,'%u PCs: [in-fold: %.1f%%] %.1f +/- %.1f%%\n',...
                        i,inSampleStats(i,1),cfnRatePCs(i),stdAcc(i));
    else
        [cfnRatePCs(i),stdAcc(i)] = ComputeCVAccuracies(topPCs,TimeSeries.Group,cfnParams);
        fprintf(1,'%u PCs: %.1f +/- %.1f%%\n',i,cfnRatePCs(i),stdAcc(i));
    end
end

%-------------------------------------------------------------------------------
%% Compute accuracy in low-dimensional space
% (first time hitting max in-sample accuracy if doComputeAtMaxInSample is true)
%-------------------------------------------------------------------------------
if doComputeAtMaxInSample
    if any(inSampleStats(:,1)==100)
        theNumPCs = find(inSampleStats(:,1)==100,1,'first');
        accAtPerfectInSample = cfnRatePCs(theNumPCs);
        stdAtPerfectInSample = stdAcc(theNumPCs);
        allFoldsAtPerfectInSample = allFoldPCs{theNumPCs};
    else
        % Keep searching
        fprintf(1,'Let''s keep searching for perfect in-sample accuracy...\n');
        PC_inSampleMeanAcc = 0;
        i = numPCs;
        while PC_inSampleMeanAcc < 100
            i = i + 1;
            topPCs = pcScore(:,1:i);
            [PCacc_mean,PCAcc_std,PC_inSampleMeanAcc,PC_allFoldData] = ComputeCVAccuracies(topPCs,TimeSeries.Group,cfnParams,true);
            fprintf(1,'%u PCs: [%.1f%%] %.1f +/- %.1f%%\n',i,PC_inSampleMeanAcc(1),PCacc_mean,PCAcc_std);
        end
        theNumPCs = i;
        accAtPerfectInSample = PCacc_mean;
        stdAtPerfectInSample = PCAcc_std;
        allFoldsAtPerfectInSample = PC_allFoldData;
    end
else
    accAtPerfectInSample = NaN;
    stdAtPerfectInSample = NaN;
    allFoldsAtPerfectInSample = NaN;
end

%-------------------------------------------------------------------------------
%% Comparison to all features
%-------------------------------------------------------------------------------
fprintf('Now with all %u features for comparison...\n',numFeatures)
if cfnParams.computePerFold
    [cfnRateAll,stdAll,inSampleAll,allFoldData] = ComputeCVAccuracies(TS_DataMat,TimeSeries.Group,cfnParams,true);
else
    [cfnRateAll,stdAll] = ComputeCVAccuracies(TS_DataMat,TimeSeries.Group,cfnParams);
end

%-------------------------------------------------------------------------------
%% Simple null estimate
%-------------------------------------------------------------------------------
numNulls = 1000;
nullStatsAll = GiveMeSimpleNullStats(TimeSeries.Group,numNulls,cfnParams,false);
nullStatsFoldMeans = mean(nullStatsAll,2);

%-------------------------------------------------------------------------------
%% Plot
%-------------------------------------------------------------------------------
plotColors = {[0,114,178]/255,[213,94,0]/255,[0,158,115]/255,[204,121,167]/255,...
                [230,159,0]/255,[86,180,233]/255,[240,215,66]/255,[0,0,0]/255};
lineWidth = 2;
f = figure('color','w');
f.Position(3:4) = [506,324];
ax = subplot(1,4,1:2);
hold('on')

%-------------------------------------------------------------------------------
% All-feature out-of-sample:
l_allfeats = yline(cfnRateAll,'-','color',plotColors{4},'LineWidth',lineWidth,...
            'Label',sprintf('All %u features (%.1f +/- %.1f%%)',numFeatures,cfnRateAll,stdAll));
legendEntryAllFeatures = sprintf('All %u features (%.1f%%)',numFeatures,cfnRateAll);

% All-feature in-sample:
if cfnParams.computePerFold
    legendEntryAllFeatures_inSample = sprintf('All features (in-fold): %.1f%%',inSampleAll(1));
    yline(inSampleAll(1),'--','color',brighten(plotColors{3},-0.5),'LineWidth',lineWidth,...
                    'Label',legendEntryAllFeatures_inSample)

end

%-------------------------------------------------------------------------------
% Naive null estimate
l_naive_null = yline(mean(nullStatsFoldMeans),'-','color',ones(1,3)*0.5,'LineWidth',lineWidth,...
                        'Label','Naive null');
l_naive_null_over = yline(quantile(nullStatsFoldMeans,0.95),'--','color',ones(1,3)*0.5,'LineWidth',lineWidth,...
                        'Label','Naive null 95% quantile');

%-------------------------------------------------------------------------------
% Out-of-fold accuracy at perfect in-fold accuracy:
if doComputeAtMaxInSample
    yline(accAtPerfectInSample,'-','color','r','LineWidth',lineWidth,...
            'Label',sprintf('At perfect in-fold acc (%.1f +/- %.1f%%)',...
                        accAtPerfectInSample,stdAtPerfectInSample));
    if theNumPCs < numPCs
        xline(theNumPCs,'--r');
    end
end

% Leading PCA accuracies:
l_testSet = plot(1:numPCs,cfnRatePCs,'o-k','LineWidth',lineWidth);
legendEntryPCs = sprintf('PCs (%.1f-%.1f%%)',min(cfnRatePCs),max(cfnRatePCs));
if cfnParams.computePerFold
    l_trainSet = plot(1:numPCs,inSampleStats(:,1),'o-','Color',plotColors{2},'LineWidth',lineWidth);
    legendEntryPCs_inSample = sprintf('In-sample (%.1f-%.1f%%)',...
                                min(inSampleStats(:,1)),max(inSampleStats(:,1)));
end

% Error bars:
plot(1:numPCs,cfnRatePCs + stdAcc,'-','Color',0.5*ones(1,3),'LineWidth',lineWidth/2);
plot(1:numPCs,cfnRatePCs - stdAcc,'-','Color',0.5*ones(1,3),'LineWidth',lineWidth/2);

% All-feature accuracies:
% yline(cfnRateAll - stdAll,':','color',plotColors{4},'LineWidth',lineWidth)
% yline(cfnRateAll + stdAll,':','color',plotColors{4},'LineWidth',lineWidth)

legend([l_testSet,l_trainSet],{legendEntryPCs,legendEntryPCs_inSample},...
        'Location','SouthEast')

ax.XTick = 1:numPCs;
xlabel('Number of PCs');
ylabel('Classification accuracy (%)')

titleText = sprintf('Classification rate (%u-class) using %u-fold %s',...
            cfnParams.numClasses,cfnParams.numFolds,cfnParams.classifierText);
title(titleText,'interpreter','none')

%-------------------------------------------------------------------------------
% Violin plots of individual fold accuracies
%-------------------------------------------------------------------------------
title(sprintf('All (%u) individual folds',cfnParams.numRepeats*cfnParams.numFolds))
ax_bar_folds = subplot(1,4,4);

if doComputeAtMaxInSample
    allFeatures_allTrainFolds = squeeze(allFoldData(:,1,:));
    allFeatures_allTestFolds = squeeze(allFoldData(:,2,:));
    allFoldsAtPerfectInSample_allTrainFolds = squeeze(allFoldsAtPerfectInSample(:,1,:));
    allFoldsAtPerfectInSample_allTestFolds = squeeze(allFoldsAtPerfectInSample(:,2,:));
    nullStatsAllFolds = nullStatsAll(:);
    accData = {nullStatsAllFolds,allFeatures_allTrainFolds(:),allFeatures_allTestFolds(:),...
            allFoldsAtPerfectInSample_allTrainFolds(:),allFoldsAtPerfectInSample_allTestFolds(:)};
    extraParams = struct();
    plotColors{1} = ones(1,3)*0.5;
    extraParams.theColors = plotColors;
    BF_ViolinPlot(accData,[],[],false,extraParams);

    ax_bar_folds.XTick = [1:5+0.5];
    xTickLabels = {'Naive-shuffle null',sprintf('all %u features train folds',numFeatures),...
                    'all features test folds',...
                    sprintf('perfect-in-fold (%u PCs) train folds',theNumPCs),...
                    sprintf('perfect-in-fold (%u PCs) test folds',theNumPCs)};

    % Bar + error bars:
    % hold('on')
    % b = bar(accData,'FaceColor','flat');
    % b.CData(1,:) = plotColors{4};
    % b.CData(2,:) = brighten(plotColors{3},-0.5);
    % b.CData(3,:) = [1,0,0];
    % stdData = [stdAll,0,stdAtPerfectInSample];
    % er = errorbar(1:3,accData,stdData);
    % er.Color = [0 0 0];
    % er.LineStyle = 'none';
else
    hold('on')
    b = bar([cfnRateAll,inSampleAll(1)],'FaceColor','flat');
    b.CData(1,:) = plotColors{4};
    b.CData(2,:) = brighten(plotColors{3},-0.5);
    ax_bar_folds.XTick = 1:2;
    xTickLabels = {'all','all (in-fold)'};
end
ax_bar_folds.XTickLabel = xTickLabels;
ax_bar_folds.YLim = ax.YLim;
linkaxes([ax,ax_bar_folds],'y')

%-------------------------------------------------------------------------------
% Violin plots of fold-average accuracies
%-------------------------------------------------------------------------------
title(sprintf('%u repeats of %u-fold averages',cfnParams.numRepeats,cfnParams.numFolds))
ax_bar = subplot(1,4,3);
hold('on')
allFeatures_allTrainFoldMeans = mean(squeeze(allFoldData(:,1,:)),1);
allFeatures_allTestFoldMeans = mean(squeeze(allFoldData(:,2,:)),1);
if ~isnan(allFoldsAtPerfectInSample)
    allFoldsAtPerfectInSample_allTrainFoldMeans = mean(squeeze(allFoldsAtPerfectInSample(:,1,:)),1);
    allFoldsAtPerfectInSample_allTestFoldMeans = mean(squeeze(allFoldsAtPerfectInSample(:,2,:)),1);
    accData = {nullStatsFoldMeans,allFeatures_allTrainFoldMeans,allFeatures_allTestFoldMeans,...
        allFoldsAtPerfectInSample_allTrainFoldMeans,allFoldsAtPerfectInSample_allTestFoldMeans};
else
    accData = {nullStatsFoldMeans,allFeatures_allTrainFoldMeans,allFeatures_allTestFoldMeans};
end

extraParams = struct();
plotColors{1} = ones(1,3)*0.5;
extraParams.theColors = plotColors;
BF_ViolinPlot(accData,[],[],false,extraParams);
l_naive_null_over = yline(quantile(nullStatsFoldMeans,0.95),'--','color',ones(1,3)*0.5,'LineWidth',lineWidth,...
                        'Label','Naive null 95% quantile');

if ~isnan(allFoldsAtPerfectInSample)
    ax_bar.XTick = (1:5+0.5);
    xTickLabels = {'Naive-shuffle null',sprintf('all %u features train-fold-means',numFeatures),...
                'all features test-fold-means',...
                sprintf('perfect-in-fold (%u PCs) train-fold-means',theNumPCs),...
                sprintf('perfect-in-fold (%u PCs) test-fold-means',theNumPCs)};
else
    ax_bar.XTick = (1:3+0.5);
    xTickLabels = {'Naive-shuffle null',sprintf('all %u features train-fold-means',numFeatures),...
                'all features test-fold-means'};
end
ax_bar.XTickLabel = xTickLabels;
ax_bar.YLim = ax.YLim;
linkaxes([ax_bar,ax_bar_folds],'y')

end
