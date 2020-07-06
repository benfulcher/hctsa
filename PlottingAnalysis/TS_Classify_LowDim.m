function TS_Classify_LowDim(whatData,cfnParams,numPCs)
% TS_Classify_LowDim compare performance of reduced PCs from the data matrix
%-------------------------------------------------------------------------------
% 'numPCs', investigate classification using up to this many PCs of the data
%              matrix (default: 0).

if nargin < 3
% numPCs:
    default_numPCs = 5;
    check_numPCs = @(x) isnumeric(x);
    addParameter(inputP,'numPCs',default_numPCs,check_numPCs);
end

%-------------------------------------------------------------------------------
% Check for NaNs in data matrix
%-------------------------------------------------------------------------------
if any(isnan(TS_DataMat(:)))
    error(['Cannot compute PCs of data matrix containing NaNs...\n' ...
            '(to compute PCs, re-run TS_Normalize to filter out all NaNs)'])
end

%-------------------------------------------------------------------------------
% Compute top X PCs of the data matrix:
%-------------------------------------------------------------------------------
fprintf('Computing top %u PCs...',numPCs)
[~,pcScore,~,~,~] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
fprintf(' Done.\n')
numPCs = min(numPCs,size(pcScore,2)); % sometimes lower than attempted 10

%-------------------------------------------------------------------------------
% Compute cumulative performance in PC space:
%-------------------------------------------------------------------------------
cfnRate = zeros(numPCs,1);
fprintf('Computing classification rates keeping top 1-%u PCs...\n',numPCs)
for i = 1:numPCs
    cfnRate(i) = GiveMeFoldLosses(pcScore(:,1:i),TimeSeries.Group);
    fprintf(1,'%u PCs:  %.3f%%\n',i,cfnRate(i));
    % fprintf(1,'%u PCs:   %.3f +/- %.3f%%\n',i,cfnRate(i,1),cfnRate(i,2));
end

%-------------------------------------------------------------------------------
% Plot
%-------------------------------------------------------------------------------
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
            cfnParams.numClasses,cfnParams.numFolds,cfnParams.whatClassifier);
title(titleText,'interpreter','none')
