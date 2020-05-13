function fn_handle = GiveMeFunctionHandle(whatClassifier,numClasses,whatLoss,reWeight)
% GiveMeFunctionHandle    Returns a function handle for classification
% accuracy and model
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
    reWeight = 1;
end
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Set the function handle to compute the accuracy/loss measure:
%-------------------------------------------------------------------------------
if strcmp(whatClassifier,'fast_linear')
  
    fn_loss = @(yTest,yPredict) BF_lossFunction(yTest,yPredict,whatLoss,numClasses);
    
    % Need to force the prior to uniform otherwise we get poor results
    % (makes it equiv. to 'classify' function)
    fn_handle_mdl = @(XTrain, yTrain) fitcdiscr(XTrain,yTrain,'Prior','uniform');
    
    fn_handle_ac = @(XTest,yTest,Mdl) fn_loss(yTest,predict(Mdl,XTest));
    fn_handle = @(XTrain,yTrain,XTest,yTest) fcnHandle(XTrain,yTrain,XTest,yTest,fn_handle_mdl,fn_handle_ac,whatLoss);
    return
end

%-------------------------------------------------------------------------------
% Binary model (easier), we can do it in one line:
%-------------------------------------------------------------------------------
if numClasses==2

    switch whatClassifier
    case 'knn'
        fn_handle_mdl = @(XTrain,yTrain,XTest,yTest) fitcknn(XTrain,yTrain);
    case 'tree'
        fn_handle_mdl = @(XTrain,yTrain,XTest,yTest) fitctree(XTrain,yTrain);
    case {'linear','linclass'}
        if reWeight
            fn_handle_mdl = @(XTrain,yTrain) fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on','Weights',InverseProbWeight(yTrain));
        else
            fn_handle_mdl = @(XTrain,yTrain) fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on');
        end
    case 'svm_linear'
        if reWeight
            fn_handle_mdl = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',InverseProbWeight(yTrain));
        else
            fn_handle_mdl = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','linear');
        end
    case 'svm_rbf'
        if reWeight
            fn_handle_mdl = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',InverseProbWeight(yTrain));
        else
            fn_handle_mdl = @(XTrain,yTrain) fitcsvm(XTrain,yTrain,'KernelFunction','rbf');
        end
    case 'diaglinear'
        if reWeight
            fn_handle_mdl = @(XTrain,yTrain) fitcnb(XTrain,yTrain,'Weights',InverseProbWeight(yTrain));
        else
            fn_handle_mdl = @(XTrain,yTrain) fitcnb(XTrain,yTrain);
        end
    otherwise
        error('Unknown classifier label: ''%s''',whatClassifier);
    end
    
    % Set the loss function:
    fn_loss = @(yTest,yPredict) BF_lossFunction(yTest,yPredict,whatLoss,numClasses);
    fn_handle_ac = @(XTest,yTest,Mdl) fn_loss(yTest,predict(Mdl,XTest));
    
    fn_handle = @(XTrain,yTrain,XTest,yTest) fcnHandle(XTrain,yTrain,XTest,yTest,fn_handle_mdl,fn_handle_ac,whatLoss);
else
    % Multiclass, have to use the Cfn function (which is slower to use as a function handle)
    fn_handle = @(XTrain,yTrain,XTest,yTest) GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,0,whatLoss,reWeight);
end

  function [accuracy,Mdl,whatLoss] = fcnHandle(XTrain, yTrain, xTest, yTest, fn_handle_mdl, fn_handle_ac, whatLoss)
    Mdl = fn_handle_mdl(XTrain, yTrain);
    accuracy = fn_handle_ac(xTest,yTest,Mdl);
  end

end
