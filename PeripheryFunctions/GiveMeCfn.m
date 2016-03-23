function [accuracy,balancedAccuracy,Mdl] = GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,beVerbose,doSumLoss)
% GiveMeFunctionHandle    Returns classification results from training a classifier on training/test data
%
%---INPUTS:
% whatClassifier -- a type of classifier to use (default: 'fast_linear')
% XTrain -- training data matrix
% yTrain -- training data labels
% XTest -- testing data matrix
% yTest -- testing data labels
% numClasses -- number of classes of labels
% beVerbose -- whether to give text output of progress
% doSumLoss -- whether to return a function measuring a total loss measure (1),
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
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 6
    numClasses = length(unique(yTrain));
end
if nargin < 7
    beVerbose = 0;
end
if nargin < 8
    doSumLoss = 0;
end

%-------------------------------------------------------------------------------
% Set the classification model:
%-------------------------------------------------------------------------------
if numClasses==2
    % Binary model (easier):
    switch whatClassifier
    case 'knn'
        Mdl = fitcknn(XTrain,yTrain,'NumNeighbors',3);
    case 'tree'
        Mdl = fitctree(XTrain,yTrain)
    case {'linear','linclass'}
        Mdl = fitcdiscr(XTrain,yTrain);
    case 'svm_linear'
        % Weight observations by inverse class probability:
        Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',ObsWeight(yTrain));
    case 'svm_rbf'
        % Weight observations by inverse class probability:
        Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',ObsWeight(yTrain));
    case 'diaglinear'
        Mdl = fitcnb(XTrain,yTrain);
    end
else
    % (Multiclass model; harder)
    switch whatClassifier
    case 'knn'
       t = templateKNN('NumNeighbors',3,'Standardize',1);
       if beVerbose, fprintf(1,'Using knn classifier\n'); end
    case 'tree'
       t = templateTree('Surrogate','off');
       if beVerbose, fprintf(1,'Using a tree classifier\n'); end
    case {'linear','linclass'}
        t = templateDiscriminant('DiscrimType','linear');
        if beVerbose, fprintf(1,'Using a linear classifier\n'); end
    case 'diaglinear'
        % Naive Bayes classifier:
        t = templateNaiveBayes('DistributionNames','normal');
        if beVerbose, fprintf(1,'Using a naive Bayes classifier\n'); end
    case 'svm_linear'
        t = templateSVM('KernelFunction','linear');
        if beVerbose, fprintf(1,'Using a linear svm classifier\n'); end
    case 'svm_rbf'
        t = templateSVM('KernelFunction','rbf');
        if beVerbose, fprintf(1,'Using a rbf svm classifier\n'); end
    end

    %-------------------------------------------------------------------------------
    % Fit the model:
    if ismember(whatClassifier,{'svm_linear','svm_rbf'})
        % Weight across potential class imbalance:
        Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'Weights',ObsWeight(yTrain));
    else
        Mdl = fitcecoc(XTrain,yTrain,'Learners',t); %,'KFold',10);
    end
end

%-------------------------------------------------------------------------------
% Evaluate performance on test data:
%-------------------------------------------------------------------------------
yPredict = predict(Mdl,XTest);
if doSumLoss
    accuracy = sum(yTest ~= yPredict); % sum of errors
else
    accuracy = mean(yTest == yPredict); % overall classification rate
end
if nargout > 1
    balancedAccuracy = mean(arrayfun(@(x) mean(yPredict(yTest==x)~=yTest(yTest==x)),unique(yTest)));
end

% switch whatClassifier
% case {'linear','linclass'}
%     fprintf(1,'A linear classifier\n');
%     Classify_fn_label = @(XTrain,yTrain,Xtest)(classify(Xtest,XTrain,yTrain,'linear'));
%     Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
%                     sum(ytest ~= classify(Xtest,XTrain,yTrain,'linear'));
% case 'diaglinear' % Naive Bayes
%     fprintf(1,'A Naive Bayes classifier\n');
%     Classify_fn_label = @(XTrain,yTrain,Xtest)(classify(Xtest,XTrain,yTrain,'diaglinear'));
%     Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
%                     sum(ytest ~= classify(Xtest,XTrain,yTrain,'diaglinear'));
% case {'svm','svmlinear'}
%     fprintf(1,'A linear support vector machine\n');
%     Classify_fn_label = @(XTrain,yTrain,Xtest) ...
%                     svmclassify(svmtrain(XTrain,yTrain, ...
%                                 'Kernel_Function','linear'),Xtest);
%     Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
%                     sum(ytest ~= svmclassify(svmtrain(XTrain,yTrain, ...
%                                 'Kernel_Function','linear'),Xtest));
% otherwise
%     error('Unknown classification method ''%s''',criterion)
% end

%-------------------------------------------------------------------------------
function weights = ObsWeight(Labels)

	classNames = unique(Labels);
	numClasses = length(classNames);

	weights = zeros(size(Labels));

	for x = 1:numClasses
		isClass = (Labels == classNames(x));
		weights(isClass) = length(Labels)/sum(isClass);
	end
end


end
