function [accuracy,accuracyStd,CVMdl,inSampleAccStd] = ComputeCVAccuracies(dataMatrix,dataLabels,cfnParams)
% Summarizes the output of GiveMeCfn as mean/std in folds, and works with single inputs to GiveMeCfn
% (for the case of CV where the training/test split is specified through folds)

[foldAcc,CVMdl] = GiveMeCfn(dataMatrix,dataLabels,dataMatrix,dataLabels,cfnParams);

if cfnParams.computePerFold
    accuracy = mean(foldAcc{2}); % test folds
    accuracyStd = std(foldAcc{2}); % test folds
    inSampleAccStd = [mean(foldAcc{1}),std(foldAcc{1})]; % training folds
else
    accuracy = foldAcc;
    accuracyStd = NaN;
end

end
