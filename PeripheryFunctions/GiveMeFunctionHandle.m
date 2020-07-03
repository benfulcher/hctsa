function fn_handle = GiveMeFunctionHandle(whatClassifier,numClasses,whatLoss,doReweight)
% GiveMeFunctionHandle  A function handle for classification accuracy and model
%
%---INPUTS:
% whatClassifier -- a type of classifier to use (default: 'fast_linear')
% numClasses -- the number of classes
% whatLoss -- a custom loss function:
%               (*) 'acc': classification rate (%)
%               (*) 'balancedAcc': mean classification rate of each class (%,
%                                        same as acc when balanced classes)
%               (*) 'sumLoss': total number of misclassifications
% doReweight -- whether to reweight observations for compatible classifiers (for
%               class imbalanced problems)

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
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
    doReweight = true;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Aim is to set the function handle for computing the accuracy/loss measure
% as a function of: @(XTrain,yTrain,XTest,yTest)
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% This is quick for the 'fast_linear' classifier:
%-------------------------------------------------------------------------------
if strcmp(whatClassifier,'fast_linear')

    % Function for fitting the model
    % (forcing a uniform prior makes this equivalent to the 'classify' function)
    fn_FitModel = @(XTrain,yTrain) fitcdiscr(XTrain,yTrain,'Prior','uniform');

    % Function for computing the loss:
    fn_Loss = @(yTest,yPredict) BF_LossFunction(yTest,yPredict,whatLoss,numClasses);

    % Combined function for computing loss and model:
    fn_handle = @(XTrain,yTrain,XTest,yTest) makeClassifyFun(XTrain,yTrain,XTest,yTest,fn_FitModel,fn_Loss,whatLoss);
    return
end

%-------------------------------------------------------------------------------
% Other classifiers:
%-------------------------------------------------------------------------------
if numClasses==2
    % Binary model is easier, we can do it in one line:

    switch whatClassifier
    case 'knn'
        fn_FitModel = @(XTrain,yTrain,XTest,yTest) fitcknn(XTrain,yTrain);
    case 'tree'
        fn_FitModel = @(XTrain,yTrain,XTest,yTest) fitctree(XTrain,yTrain);
    case {'linear','linclass'}
        if doReweight
            fn_FitModel = @(XTrain,yTrain) fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on','Weights',InverseProbWeight(yTrain));
        else
            fn_FitModel = @(XTrain,yTrain) fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on');
        end
    case 'svm_linear'
        if doReweight
            fn_FitModel = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','linear',...
                                                    'Weights',InverseProbWeight(yTrain));
        else
            fn_FitModel = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','linear');
        end
    case 'svm_rbf'
        if doReweight
            fn_FitModel = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',InverseProbWeight(yTrain));
        else
            fn_FitModel = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','rbf');
        end
    case 'diaglinear'
        if doReweight
            fn_FitModel = @(XTrain,yTrain) fitcnb(XTrain,yTrain,'Weights',InverseProbWeight(yTrain));
        else
            fn_FitModel = @(XTrain,yTrain) fitcnb(XTrain,yTrain);
        end
    otherwise
        error('Unknown classifier: ''%s''',whatClassifier);
    end

    % Set the loss function:
    fn_Loss = @(yTest,yPredict) BF_LossFunction(yTest,yPredict,whatLoss,numClasses);

    % Combine into a single function handle:
    fn_handle = @(XTrain,yTrain,XTest,yTest) makeClassifyFun(XTrain,yTrain,XTest,yTest,fn_FitModel,fn_Loss,whatLoss);
else
    % Multiclass, have to use the Cfn function (which is slower to use as a function handle)
    fn_handle = @(XTrain,yTrain,XTest,yTest) GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,0,whatLoss,doReweight);
end

%-------------------------------------------------------------------------------
function [accuracy,Mdl,whatLoss] = makeClassifyFun(XTrain,yTrain,xTest,yTest,fn_FitModel,fn_Loss,whatLoss)
    % A function that:
    % 1. trains the model (Mdl)
    % 2. computes its accuracy on test data
    % 3. returns the type of loss function used
    Mdl = fn_FitModel(XTrain,yTrain);
    accuracy = fn_Loss(yTest,predict(Mdl,XTest));
end

end
