function [accuracyMean,accuracyStd,inSampleAccStd,allFoldData,CVMdl] = ComputeCVAccuracies(dataMatrix,dataLabels,cfnParams,doSummarize)
% ComputeCVAccuracies
% Summarizes the output of GiveMeCfn as mean/std in folds, and works with single inputs to GiveMeCfn
% (for the case of CV where the training/test split is specified through folds)

%-------------------------------------------------------------------------------
%% Check inputs and set defaults
%-------------------------------------------------------------------------------
if nargin < 4
    doSummarize = true;
end

%-------------------------------------------------------------------------------
%% Fit the CV model
%-------------------------------------------------------------------------------
% Initialize
CVMdl = cell(cfnParams.numRepeats,1);
accuracyMean = zeros(cfnParams.numRepeats,1);
accuracyStd = zeros(cfnParams.numRepeats,1);
if cfnParams.computePerFold
    inSampleAccStd = zeros(cfnParams.numRepeats,2);
else
    % Not appropriate: set to NaN in case output specified
    inSampleAccStd = NaN;
end
allFoldData = zeros(cfnParams.numFolds,2,cfnParams.numRepeats); % folds, train/test, repeats

% Fit and evaluate the model using a CV partition numRepeats times:
for r = 1:cfnParams.numRepeats
    [foldAcc,CVMdl{r}] = GiveMeCfn(dataMatrix,dataLabels,dataMatrix,dataLabels,cfnParams);

    if cfnParams.computePerFold
        allFoldData(:,1,r) = foldAcc{1}; % training folds
        allFoldData(:,2,r) = foldAcc{2}; % test folds
        % Summaries as mean/std over folds:
        accuracyMean(r) = mean(foldAcc{2}); % test folds
        accuracyStd(r) = std(foldAcc{2}); % test folds
        inSampleAccStd(r,:) = [mean(foldAcc{1}),std(foldAcc{1})]; % training folds
    else
        % Already outputs a mean over folds:
        accuracyMean(r) = foldAcc;
        accuracyStd(r) = NaN;
    end
end

%-------------------------------------------------------------------------------
%% Summarize quantities as the mean across repeats
%-------------------------------------------------------------------------------
if doSummarize
    accuracyMean = mean(accuracyMean,1);
    accuracyStd = mean(accuracyStd,1);
    inSampleAccStd = mean(inSampleAccStd,1);
    % allFoldData = mean(allFoldData,3);
end

end
