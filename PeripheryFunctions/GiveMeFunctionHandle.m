function fn_handle = GiveMeFunctionHandle(cfnParams)
% GiveMeFunctionHandle    Returns a function handle for computing classification
%                               accuracy
%
%---INPUTS:
% cfnParams, parameters for the classification, e.g., set with
%           GiveMeDefaultClassificationParams.

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
% Set the function handle to compute the accuracy/loss measure:
%-------------------------------------------------------------------------------
if strcmp(cfnParams.whatClassifier,'fast_linear')
    fn_loss = @(yTest,yPredict) BF_LossFunction(yTest,yPredict,cfnParams.whatLoss,cfnParams.classLabels);
    fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,classify(XTest,XTrain,yTrain,'linear'));
    return
end

%-------------------------------------------------------------------------------
% Binary model (easier), we can do it in one line:
%-------------------------------------------------------------------------------
if cfnParams.numClasses==2
    % Set the loss function:
    fn_loss = @(yTest,yPredict) BF_LossFunction(yTest,yPredict,cfnParams.whatLoss,cfnParams.classLabels);

    switch cfnParams.whatClassifier
    case 'knn'
        fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcknn(XTrain,yTrain),XTest));
    case 'tree'
        fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitctree(XTrain,yTrain),XTest));
    case {'linear','linclass'}
        if cfnParams.doReweight
            fn_handle = @(XTrain,yTrain,XTest,yTest) ...
                        fn_loss(yTest,predict(fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) ...
                        fn_loss(yTest,predict(fitcdiscr(XTrain,yTrain,'FillCoeffs','off',...
                                    'SaveMemory','on'),XTest));
        end
    case 'svm_linear'
        if cfnParams.doReweight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','linear','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','linear'),XTest));
        end
    case 'svm_rbf'
        if cfnParams.doReweight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','rbf','Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcsvm(XTrain,yTrain,...
                                    'KernelFunction','rbf'),XTest));
        end
    case 'diaglinear'
        if cfnParams.doReweight
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcnb(XTrain,yTrain,...
                                                    'Weights',InverseProbWeight(yTrain)),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) fn_loss(yTest,predict(fitcnb(XTrain,yTrain),XTest));
        end
    otherwise
        error('Unknown classifier: ''%s''',cfnParams.whatClassifier);
    end
else
    % Multiclass, have to use the Cfn function (which is slower to use as a function handle)
    fn_handle = @(XTrain,yTrain,XTest,yTest) ...
        GiveMeCfn(XTrain,yTrain,XTest,yTest,cfnParams,false);
end

end
