function [accuracyMean,accuracyStd,inSampleAccStd,allFoldData,CVMdl] = ComputeCVAccuracies(dataMatrix,dataLabels,cfnParams,doSummarize)
% ComputeCVAccuracies   Classification performance from repeated cross-validation.
%
% Fits the classifier specified in cfnParams (using GiveMeCfn) cfnParams.numRepeats
% times, and summarizes its performance (cfnParams.whatLoss) across folds.
%
%---INPUTS:
% dataMatrix, the data matrix (observations x features)
% dataLabels, the class label of each observation
% cfnParams, classification parameters (cf. GiveMeDefaultClassificationParams)
% doSummarize, [default: true] whether to average the outputs across repeats
%
%---OUTPUTS:
% (one row per repeat, or the mean across repeats if doSummarize = true)
% accuracyMean, the performance on the test folds: the mean across folds if
%               cfnParams.computePerFold, otherwise computed once from the
%               predictions of all test folds. If cfnParams.numFolds = 0 (no
%               cross-validation), the in-sample performance.
% accuracyStd, the standard deviation of the performance across test folds
%               (NaN unless cfnParams.computePerFold and cfnParams.numFolds > 0).
% inSampleAccStd, [mean,std] of the performance across training folds (NaN
%               unless cfnParams.computePerFold and cfnParams.numFolds > 0).
% allFoldData, the performance in each fold, as numFolds x 2 (training, test) x
%               numRepeats (NaN unless cfnParams.computePerFold; empty if
%               cfnParams.numFolds = 0). It is never averaged across repeats.
% CVMdl, the fitted model from each repeat (cell): cross-validated models if
%               cfnParams.numFolds > 0, otherwise models trained on all data.

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 4
    doSummarize = true;
end

% Results are only computed per fold when there are folds:
perFold = cfnParams.computePerFold && (cfnParams.numFolds > 0);
cfnParams.computePerFold = perFold;

%-------------------------------------------------------------------------------
%% Fit the CV model
%-------------------------------------------------------------------------------
% Initialize
CVMdl = cell(cfnParams.numRepeats,1);
accuracyMean = zeros(cfnParams.numRepeats,1);
accuracyStd = nan(cfnParams.numRepeats,1);
inSampleAccStd = nan(cfnParams.numRepeats,2);
allFoldData = nan(cfnParams.numFolds,2,cfnParams.numRepeats); % folds, train/test, repeats

% Fit and evaluate the model using a CV partition numRepeats times:
for r = 1:cfnParams.numRepeats
    [foldAcc,CVMdl{r}] = GiveMeCfn(dataMatrix,dataLabels,dataMatrix,dataLabels,cfnParams);

    if perFold
        allFoldData(:,1,r) = foldAcc{1}; % training folds
        allFoldData(:,2,r) = foldAcc{2}; % test folds
        % Summaries as mean/std over folds:
        accuracyMean(r) = mean(foldAcc{2}); % test folds
        accuracyStd(r) = std(foldAcc{2}); % test folds
        inSampleAccStd(r,:) = [mean(foldAcc{1}),std(foldAcc{1})]; % training folds
    else
        % Already outputs a mean over folds (or the in-sample value without folds):
        accuracyMean(r) = foldAcc;
    end
end

%-------------------------------------------------------------------------------
%% Summarize quantities as the mean across repeats
%-------------------------------------------------------------------------------
if doSummarize
    accuracyMean = mean(accuracyMean,1);
    accuracyStd = mean(accuracyStd,1);
    inSampleAccStd = mean(inSampleAccStd,1);
end

end
