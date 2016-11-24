function fn_handle = GiveMeFunctionHandle(whatClassifier,numClasses,whatLoss,reWeight)
% GiveMeFunctionHandle    Returns a function handle for classification accuracy
%
%---INPUTS:
% whatClassifier -- a type of classifier to use (default: 'fast_linear')
% numClasses -- the number of classes
% whatLoss -- a custom loss function:
%               (*) 'acc': classification rate (%)
%               (*) 'balancedAcc': mean classification rate of each class (%,
%                                        same as acc when balanced classes)
%               (*) 'sumLoss': total number of misclassifications
% reWeight -- whether to reweight observations for compatible classifiers (for
%               class imbalanced problems)

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
if nargin < 1
    whatClassifier = 'fast_linear';
end
if nargin < 2
    warning('How many classes?! Assuming a multiclass problem')
    numClasses = 3; % assume multiclass
end
if nargin < 3
    whatLoss = 'acc'; % total accuracy (by default)
end
if nargin < 4
    reWeight = 1;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Set the function handle to compute the accuracy/loss measure:
%-------------------------------------------------------------------------------
if strcmp(whatClassifier,'fast_linear')
    fn_loss = @(yTest,yPredict) BF_lossFunction(yTest,yPredict,whatLoss,numClasses);
    fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,classify(XTest,XTrain,yTrain,'linear'));
    return
end

%-------------------------------------------------------------------------------
% Binary model (easier), we can do it in one line:
%-------------------------------------------------------------------------------
if numClasses==2
    % Set the loss function:
    fn_loss = @(yTest,yPredict) BF_lossFunction(yTest,yPredict,whatLoss,numClasses);

    switch whatClassifier
    case 'knn'
        fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcknn(XTrain,yTrain),XTest));
    case 'tree'
        fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitctree(XTrain,yTrain),XTest));
    case {'linear','linclass'}
        if reWeight
            fn_handle = @(XTrain,yTrain,XTest,yTest) ...
                        fn_loss(yTest,predict(fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) ...
                        fn_loss(yTest,predict(fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on'),XTest));
        end
    case 'svm_linear'
        if reWeight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','linear','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','linear'),XTest));
        end
    case 'svm_rbf'
        if reWeight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','rbf','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','rbf'),XTest));
        end
    case 'diaglinear'
        if reWeight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcnb(XTrain,yTrain,...
                                                    'Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcnb(XTrain,yTrain),XTest));
        end
    otherwise
        error('Unknown classifier label: ''%s''',whatClassifier);
    end
else
    % Multiclass, have to use the Cfn function (which is slower to use as a function handle)
    fn_handle = @(XTrain,yTrain,XTest,yTest) ...
        GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,0,whatLoss,reWeight);
end

end
