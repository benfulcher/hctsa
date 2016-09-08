function [accuracy,Mdl,whatLoss] = GiveMeCfn(whatClassifier,XTrain,yTrain,XTest,yTest,numClasses,beVerbose,whatLoss,reWeight,CVFolds)
% GiveMeCfn    Returns classification results from training a classifier on
%               training/test data
%
%---INPUTS:
% whatClassifier -- a type of classifier to use (default: 'fast_linear')
% XTrain -- training data matrix
% yTrain -- training data labels
% XTest -- testing data matrix
% yTest -- testing data labels
% numClasses -- number of classes of labels
% beVerbose -- whether to give text output of progress
% whatLoss -- what loss function to compute on the data

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
% Make sure y a column vector
if size(yTrain,1) < size(yTrain,2)
    yTrain = yTrain';
end
if nargin < 4
    XTest = [];
end
if nargin < 5
    yTest = [];
end
if nargin < 6
    numClasses = max(yTrain);
end
if nargin < 7
    beVerbose = 0;
end
if nargin < 8 || isempty(whatLoss)
    % See if it's a balanced problem, and set defaults accordingly
    yAll = [yTrain;yTest];
    classNumbers = arrayfun(@(x)sum(yAll==x),1:numClasses);
    isBalanced = all(classNumbers==classNumbers(1));
    if isBalanced
        whatLoss = 'acc';
        reWeight = 0;
    else
        whatLoss = 'balancedAcc';
        reWeight = 1;
        if beVerbose
            fprintf(1,'Unbalanced classes: using a balanced accuracy measure (& using reweighting)...\n');
        end
    end
end
if ~exist('reWeight','var')
    % Reweighted observations by inverse probability weight
    % (for class imbalanced problems)
    reWeight = 1;
end
if nargin < 10
    CVFolds = 0;
end

% Reinterpret input: "fast_linear"
if strcmp(whatClassifier,'fast_linear')
    if beVerbose
        fprintf(1,'Using diaglinear instead of fast_linear\n');
    end
    whatClassifier = 'diaglinear';
end

%-------------------------------------------------------------------------------
% Set the classification model:
%-------------------------------------------------------------------------------
if numClasses==2
    % Binary model (easier):
    switch whatClassifier
    case 'knn'
        if CVFolds > 0
            Mdl = fitcknn(XTrain,yTrain,'NumNeighbors',3,'KFold',CVFolds);
        else
            Mdl = fitcknn(XTrain,yTrain,'NumNeighbors',3);
        end
    case 'tree'
        if CVFolds > 0
            Mdl = fitctree(XTrain,yTrain,'KFold',CVFolds)
        else
            Mdl = fitctree(XTrain,yTrain)
        end
    case {'linear','linclass'}
        if CVFolds > 0
            Mdl = fitcdiscr(XTrain,yTrain,'FillCoeffs','off','SaveMemory','on','KFold',CVFolds);
        else
            Mdl = fitcdiscr(XTrain,yTrain,'FillCoeffs','off','SaveMemory','on');
        end
    case {'svm','svm_linear'}
        % Weight observations by inverse class probability:
        if reWeight
            if CVFolds > 0
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',InverseProbWeight(yTrain),'KFold',CVFolds);
            else
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',InverseProbWeight(yTrain));
            end
        else
            if CVFolds > 0
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','KFold',CVFolds);
            else
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear');
            end
        end
    case 'svm_rbf'
        % Weight observations by inverse class probability:
        if reWeight
            if CVFolds > 0
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',InverseProbWeight(yTrain),'KFold',CVFolds);
            else
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',InverseProbWeight(yTrain));
            end
        else
            if CVFolds > 0
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','KFold',CVFolds);
            else
                Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf');
            end
        end
    case 'diaglinear'
        if CVFolds > 0
            Mdl = fitcnb(XTrain,yTrain,'KFold',CVFolds);
        else
            Mdl = fitcnb(XTrain,yTrain);
        end
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
        t = templateDiscriminant('DiscrimType','linear','FillCoeffs','off','SaveMemory','on');
        if beVerbose, fprintf(1,'Using a linear classifier\n'); end
    case 'diaglinear'
        % Naive Bayes classifier:
        t = templateNaiveBayes('DistributionNames','normal');
        if beVerbose, fprintf(1,'Using a naive Bayes classifier\n'); end
    case {'svm','svm_linear'}
        t = templateSVM('KernelFunction','linear');
        if beVerbose, fprintf(1,'Using a linear svm classifier\n'); end
    case 'svm_rbf'
        t = templateSVM('KernelFunction','rbf');
        if beVerbose, fprintf(1,'Using a rbf svm classifier\n'); end
    end

    %-------------------------------------------------------------------------------
    % Fit the model:
    if ismember(whatClassifier,{'svm_linear','svm_rbf','linear','linclass','diaglinear'}) && reWeight
        % Weight across potential class imbalance:
        if CVFolds > 0
            Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'Weights',InverseProbWeight(yTrain),'KFold',CVFolds);
        else
            Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'Weights',InverseProbWeight(yTrain));
        end
    else
        % Don't reweight:
        if CVFolds > 0
            Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'KFold',CVFolds);
        else
            Mdl = fitcecoc(XTrain,yTrain,'Learners',t);
        end
    end
end

%-------------------------------------------------------------------------------
% Evaluate performance on test data:
%-------------------------------------------------------------------------------

% Predict for the test data:
if CVFolds == 0
    if isempty(XTest) || isempty(yTest)
        % No need to compute this if you're just after the model
        accuracy = [];
        balancedAccuracy = [];
        return
    else
        yPredict = predict(Mdl,XTest);
        accuracy = BF_lossFunction(yTest,yPredict,whatLoss,numClasses);
    end
else
    % Test data is mixed through the training data provided using k-fold cross validation
    % Output is the accuracy/loss measure for each fold, can mean it themselves if they want to
    yPredict = kfoldPredict(Mdl);
    accuracy = arrayfun(@(x) BF_lossFunction(yTrain(Mdl.Partition.test(x)),...
                                yPredict(Mdl.Partition.test(x)),whatLoss,numClasses),...
                                    1:CVFolds);
end


end
