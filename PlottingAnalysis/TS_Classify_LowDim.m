function TS_Classify_LowDim(whatData,cfnParams,numPCs)
% TS_Classify_LowDim compare performance of reduced PCs from the data matrix
%-------------------------------------------------------------------------------
% 'whatData', HCTSA data file (or structure).
% 'cfnParams', parameters of the classification to be performed (cf.
%                   GiveMeDefaultClassificationParams)
% 'numPCs', investigate classification using up to this many PCs of the data
%              matrix (default: 0).

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

if nargin < 1
    whatData = 'norm';
end
if nargin < 3
    numPCs = 5;
end

%-------------------------------------------------------------------------------
% Load data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
numFeatures = size(TS_DataMat,2);

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
[groupLabels,classLabels,groupLabelsInteger,numGroups] = ExtractGroupLabels(TimeSeries);
TellMeAboutLabeling(TimeSeries);
if nargin < 2
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end
TellMeAboutClassification(cfnParams);

%-------------------------------------------------------------------------------
% Check for NaNs in data matrix
%-------------------------------------------------------------------------------
if any(isnan(TS_DataMat(:)))
    warning(['Consider re-running TS_Normalize to filter out NaNs from the feature matrix)'])
end

%-------------------------------------------------------------------------------
% Compute top X PCs of the data matrix:
%-------------------------------------------------------------------------------
% Use pca to compute the first two principal components:
% (project data into space of PC scores, Y)
if ~any(isnan(TS_DataMat))
    fprintf('Computing top %u PCs...',numPCs)
    [pcCoeff,pcScore,~,~,percVar] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
else
    warning(sprintf(['Data matrix contains %.2g%% NaNs. Estimating covariances on remaining data...\n' ...
                '(Could take some time...)'],100*mean(isnan(TS_DataMat(:)))))
    [pcCoeff,pcScore,~,~,percVar] = pca(BF_NormalizeMatrix(TS_DataMat,'zscore'),...
                            'Rows','pairwise','NumComponents',numPCs);
    % If this fails (covariance matrix not positive definite), can try the
    % (...,'algorithm','als') option in pca... (or toolbox for probabilistic PCA)
end
fprintf(' Done.\n')
numPCs = min(numPCs,size(pcScore,2)); % sometimes lower than the number attempted

%-------------------------------------------------------------------------------
% Display some info about feature loading onto the reduced components:
%-------------------------------------------------------------------------------
numTopLoadFeat = min(numFeatures,20); % display this many features loading onto each PC
LowDimDisplayTopLoadings(numTopLoadFeat,numPCs,pcCoeff,pcScore,TS_DataMat,Operations);

%-------------------------------------------------------------------------------
% Compute cumulative performance in PC space:
%-------------------------------------------------------------------------------
cfnRatePCs = zeros(numPCs,1);
fprintf('Computing classification rates keeping top 1-%u PCs...\n',numPCs)
cfnParams.computePerFold = false;
for i = 1:numPCs
    cfnRatePCs(i) = GiveMeCfn(pcScore(:,1:i),TimeSeries.Group,[],[],cfnParams);
    fprintf(1,'%u PCs: %.3f%%\n',i,cfnRatePCs(i));
    % fprintf(1,'%u PCs:   %.3f +/- %.3f%%\n',i,cfnRate(i,1),cfnRate(i,2));
end

%-------------------------------------------------------------------------------
% Comparison to all features
%-------------------------------------------------------------------------------
fprintf('Now with all %u features for comparison...\n',numFeatures)
cfnRateAll = GiveMeCfn(TS_DataMat,TimeSeries.Group,[],[],cfnParams);

%-------------------------------------------------------------------------------
% Plot
%-------------------------------------------------------------------------------
plotColors = BF_GetColorMap('spectral',3,1);
lineWidth = 2;
f = figure('color','w');
f.Position(3:4) = [506,324];
hold('on')
ax = gca();
plot([1,numPCs],ones(2,1)*cfnRateAll,'--','color',plotColors{3},'LineWidth',lineWidth)
plot(1:numPCs,cfnRatePCs(:,1),'o-k','LineWidth',lineWidth)
legend(sprintf('All %u features (%.1f%%)',numFeatures,cfnRateAll),...
            sprintf('PCs (%.1f-%.1f%%)',min(cfnRatePCs),max(cfnRatePCs)))
% plot(1:numPCs,cfnRatePCs(:,1)+cfnRatePCs(:,2),':k')
% plot(1:numPCs,cfnRatePCs(:,1)-cfnRatePCs(:,2),':k')

ax.XTick = 1:numPCs;
xlabel('Number of PCs');
ylabel('Classification accuracy (%)')

titleText = sprintf('Classification rate (%u-class) using %u-fold %s classification',...
            cfnParams.numClasses,cfnParams.numFolds,cfnParams.whatClassifier);
title(titleText,'interpreter','none')

end
