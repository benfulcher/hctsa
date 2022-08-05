function nullStats = GiveMeSimpleNullStats(groupLabels,numNulls,cfnParams,meanOverFolds)
% Returns simple permutation-based null samples
%-------------------------------------------------------------------------------

% Check inputs and set defaults:
if nargin < 4
    meanOverFolds = true;
end

%-------------------------------------------------------------------------------
numDataSamples = length(groupLabels);

if meanOverFolds
    nullStats = zeros(numNulls,1);
else
    nullStats = zeros(numNulls,cfnParams.numFolds);
end

for i = 1:numNulls
    shuffledLabels = groupLabels(randperm(numDataSamples));

    % One repeat of 10-fold cross validation using the shuffled labels:
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
