function nullStats = GiveMeSimpleNullStats(groupLabels,numNulls,cfnParams,meanOverFolds)
% GiveMeSimpleNullStats  Simple permutation-based null samples of a loss
%
% Each null sample is the loss (cfnParams.whatLoss) of a random shuffle of the labels
% (with the true class frequencies) taken as the prediction, evaluated on the test
% sets of one stratified cross-validation partition into cfnParams.numFolds folds
% (a new partition for every null sample). It uses the global random number
% generator.
%
%---INPUTS:
% groupLabels, the (categorical) class label of each time series
% numNulls, the number of null samples to generate
% cfnParams, classification settings (needs numFolds (at least 2), whatLoss, classLabels)
% meanOverFolds, whether to average the loss over folds (default: true)
%
%---OUTPUT:
% nullStats, numNulls x 1 (meanOverFolds = true) or numNulls x numFolds
%-------------------------------------------------------------------------------

% Check inputs and set defaults:
if nargin < 4
    meanOverFolds = true;
end

%-------------------------------------------------------------------------------
if cfnParams.numFolds < 2
    error('The null requires cross-validation with at least 2 folds (cfnParams.numFolds = %u)',...
                        cfnParams.numFolds);
end
numDataSamples = length(groupLabels);

if meanOverFolds
    nullStats = zeros(numNulls,1);
else
    nullStats = zeros(numNulls,cfnParams.numFolds);
end

for i = 1:numNulls
    shuffledLabels = groupLabels(randperm(numDataSamples));

    % One repeat of cross validation using the shuffled labels:
    cvFolds = cvpartition(groupLabels,'KFold',cfnParams.numFolds,'Stratify',true);
    nullAcc_k = zeros(cfnParams.numFolds,1);
    for k = 1:cfnParams.numFolds
        yTrue = groupLabels(cvFolds.test(k));
        yPredict = shuffledLabels(cvFolds.test(k));
        nullAcc_k(k) = BF_LossFunction(yTrue,yPredict,cfnParams.whatLoss,cfnParams.classLabels);
    end

    % Agglomerate over folds:
    if meanOverFolds
        nullStats(i) = mean(nullAcc_k);
    else
        nullStats(i,:) = nullAcc_k;
    end
end

end
