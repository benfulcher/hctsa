function fn_handle = GiveMeFunctionHandle(whatClassifier,numClasses,sumLoss)
% GiveMeFunctionHandle    Returns a function handle for classification accuracy
%
%---INPUTS:
% whatClassifier -- a type of classifier to use (default: 'fast_linear')
% numClasses -- the number of classes
% sumLoss -- whether to return a function measuring a total loss measure (1),
%               or an accuracy rate (0, default)

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
    numClasses = 3; % multiclass
end
if nargin < 3
    sumLoss = 0;
end
%-------------------------------------------------------------------------------

if strcmp(whatClassifier,'fast_linear')
    if sumLoss
        fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= classify(XTest,XTrain,yTrain,'linear'));
    else
        fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == classify(XTest,XTrain,yTrain,'linear'));
    end
    return
end

if numClasses==2
    % Binary model (easier), we can do it in one line:
    switch whatClassifier
    case 'knn'
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitcknn(XTrain,yTrain),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitcknn(XTrain,yTrain),XTest));
        end
    case 'tree'
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitctree(XTrain,yTrain),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitctree(XTrain,yTrain),XTest));
        end
    case {'linear','linclass'}
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitcdiscr(XTrain,yTrain),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitcdiscr(XTrain,yTrain),XTest));
        end
    case 'svm_linear'
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitcsvm(XTrain,yTrain,'KernelFunction','linear'),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitcsvm(XTrain,yTrain,'KernelFunction','linear'),XTest));
        end
    case 'svm_rbf'
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitcsvm(XTrain,yTrain,'KernelFunction','rbf'),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitcsvm(XTrain,yTrain,'KernelFunction','rbf'),XTest));
        end
    case 'diaglinear'
        if sumLoss
            fn_handle = @(XTrain,yTrain,XTest,yTest) sum(yTest ~= predict(fitcnb(XTrain,yTrain),XTest));
        else
            fn_handle = @(XTrain,yTrain,XTest,yTest) 100*mean(yTest == predict(fitcnb(XTrain,yTrain),XTest));
        end
    otherwise
        error('Unknown classifier label: ''%s''',whatClassifier);
    end
else
    % Multiclass
    if sumLoss
        fn_handle = @(XTrain,yTrain,XTest,yTest) ...
            GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,0);
    else
        fn_handle = @(XTrain,yTrain,XTest,yTest) ...
            100*GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,0);
    end
end

end
