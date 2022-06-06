function [accuracy,Mdl] = GiveMeCfn(XTrain,yTrain,XTest,yTest,cfnParams,beVerbose)
% GiveMeCfn    Returns classification results from training a classifier on
%               training/test data
%
%---INPUTS:
% XTrain -- training data matrix
% yTrain -- training data labels
% XTest -- testing data matrix
% yTest -- testing data labels
% cfnParams -- parameters for classification
% beVerbose -- [default: false] whether to give text output of progress

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
%% Check Inputs:
%-------------------------------------------------------------------------------
% Ensure y is a column vector
if size(yTrain,1) < size(yTrain,2)
    yTrain = yTrain';
end
if nargin < 3
    XTest = [];
end
if nargin < 4
    yTest = [];
end
if nargin < 6
    beVerbose = false;
end

%-------------------------------------------------------------------------------
%% Define the classification model
%------------------------------------------------------------------------------
if strcmp(cfnParams.whatClassifier,'fast-linear')
    %--------------------------------------------------------------------------
    % Special case -- the `classify` function is faster than others
    %--------------------------------------------------------------------------
    if cfnParams.numFolds > 0
        if cfnParams.doReweight
            Mdl = fitcdiscr(XTrain,yTrain,'Prior','empirical',...
            'DiscrimType','linear','Weights',InverseProbWeight(yTrain),...
            'KFold',cfnParams.numFolds);
        else
            Mdl = fitcdiscr(XTrain,yTrain,'Prior','empirical',...
            'DiscrimType','linear','KFold',cfnParams.numFolds);
        end
    else
        if nargout == 1
            % It's faster to use classify if we only want the accuracy
            accuracy = BF_LossFunction(yTest,classify(yTest,XTrain,yTrain),...
                            cfnParams.whatLoss,cfnParams.classLabels);
            return;
        else
            Mdl = fitcdiscr(XTrain,yTrain,'Prior','empirical');
        end
    end
else
    if cfnParams.numClasses==2
        %----------------------------------------------------------------------
        % Special case: a two-class model
        %----------------------------------------------------------------------
        switch cfnParams.whatClassifier
        case 'knn'
            if beVerbose
                fprintf(1,'Using three neighbors for knn\n');
            end
            if cfnParams.numFolds > 0
                Mdl = fitcknn(XTrain,yTrain,'NumNeighbors',3,'KFold',cfnParams.numFolds);
            else
                Mdl = fitcknn(XTrain,yTrain,'NumNeighbors',3);
            end
        case 'tree'
            if cfnParams.numFolds > 0
                Mdl = fitctree(XTrain,yTrain,'KFold',cfnParams.numFolds);
            else
                Mdl = fitctree(XTrain,yTrain);
            end
        case {'linear','linclass'}
            if cfnParams.numFolds > 0
                Mdl = fitcdiscr(XTrain,yTrain,'FillCoeffs','off','SaveMemory','on','KFold',cfnParams.numFolds);
            else
                Mdl = fitcdiscr(XTrain,yTrain,'FillCoeffs','off','SaveMemory','on');
            end
        case 'logistic'
            if cfnParams.doReweight
                if cfnParams.numFolds > 0
                    Mdl = fitclinear(XTrain,yTrain,'Learner','logistic','Regularization','ridge',...
                                'Lambda','auto','Weights',InverseProbWeight(yTrain),'KFold',cfnParams.numFolds);
                else
                    Mdl = fitclinear(XTrain,yTrain,'Learner','logistic','Regularization','ridge',...
                                'Lambda','auto','Weights',InverseProbWeight(yTrain));
                end
            else
                if cfnParams.numFolds > 0
                    Mdl = fitclinear(XTrain,yTrain,'Learner','logistic','Regularization','ridge',...
                                'Lambda','auto','KFold',cfnParams.numFolds);
                else
                    Mdl = fitclinear(XTrain,yTrain,'Learner','logistic','Regularization','ridge',...
                                'Lambda','auto');
                end
            end
        case {'svm','svm-linear'}
            % Weight observations by inverse class probability:
            % lambdaDefault = 1/size(XTrain,1);
            if cfnParams.doReweight
                if cfnParams.numFolds > 0
                    Mdl = fitclinear(XTrain,yTrain,'Learner','svm','Lambda','auto','Regularization','ridge',...
                                'Weights',InverseProbWeight(yTrain),'KFold',cfnParams.numFolds);
                else
                    Mdl = fitclinear(XTrain,yTrain,'Learner','svm','Lambda','auto','Regularization','ridge',...
                                'Weights',InverseProbWeight(yTrain));
                end
            else
                if cfnParams.numFolds > 0
                    Mdl = fitclinear(XTrain,yTrain,'Learner','svm','Lambda','auto',...
                                        'Regularization','ridge','KFold',cfnParams.numFolds);
                else
                    Mdl = fitclinear(XTrain,yTrain,'Learner','svm','Lambda','auto','Regularization','ridge');
                end
            end
        case 'svm-linear-lowdim'
            % For low-dim datasets (e.g., using a reduced feature set), the 'traditional' linear SVM (fitcsvm)
            % can be preferable to the above (fitclinear)
            if cfnParams.doReweight
                if cfnParams.numFolds > 0
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',InverseProbWeight(yTrain),'KFold',cfnParams.numFolds);
                else
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','Weights',InverseProbWeight(yTrain));
                end
            else
                if cfnParams.numFolds > 0
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear','KFold',cfnParams.numFolds);
                else
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','linear');
                end
            end
        case 'svm-rbf'
            % Weight observations by inverse class probability:
            if cfnParams.doReweight
                if cfnParams.numFolds > 0
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',...
                                    InverseProbWeight(yTrain),'KFold',cfnParams.numFolds);
                else
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','Weights',...
                                                    InverseProbWeight(yTrain));
                end
            else
                if cfnParams.numFolds > 0
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf','KFold',cfnParams.numFolds);
                else
                    Mdl = fitcsvm(XTrain,yTrain,'KernelFunction','rbf');
                end
            end
        case 'diaglinear'
            if cfnParams.numFolds > 0
                Mdl = fitcnb(XTrain,yTrain,'KFold',cfnParams.numFolds);
            else
                Mdl = fitcnb(XTrain,yTrain);
            end
        otherwise
            error('Unknown classifier: ''%s''',cfnParams.whatClassifier);
        end
    else
        %----------------------------------------------------------------------
        % Define and fit a classification model involving 3 or more classes
        %----------------------------------------------------------------------
        % Define the model:
        switch cfnParams.whatClassifier
        case 'knn'
            t = templateKNN('NumNeighbors',3,'Standardize',true);
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
        case {'svm','svm-linear'}
            t = templateSVM('KernelFunction','linear');
            if beVerbose, fprintf(1,'Using a linear svm classifier\n'); end
        case 'svm-rbf'
            t = templateSVM('KernelFunction','rbf');
            if beVerbose, fprintf(1,'Using a rbf svm classifier\n'); end
        case 'logistic'
            error('Logistic regression is not applicable to problems with >2 classes')
        otherwise
            error('Unknown classifier: ''%s''',cfnParams.whatClassifier);
        end

        % Fit the model:
        if ismember(cfnParams.whatClassifier,{'svm-linear','svm-rbf','linear','linclass','diaglinear'}) && cfnParams.doReweight
            % Reweight to give equal weight to each class (in case of class imbalance)
            if cfnParams.numFolds > 0
                Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'Weights',InverseProbWeight(yTrain),'KFold',cfnParams.numFolds);
            else
                Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'Weights',InverseProbWeight(yTrain));
            end
        else
            % Don't reweight:
            if cfnParams.numFolds > 0
                Mdl = fitcecoc(XTrain,yTrain,'Learners',t,'KFold',cfnParams.numFolds);
            else
                Mdl = fitcecoc(XTrain,yTrain,'Learners',t);
            end
        end
    end
end

%-------------------------------------------------------------------------------
%% Evaluate performance on test data (in-sample, or from CV folds)
%-------------------------------------------------------------------------------
% Predict the test data:
if cfnParams.numFolds == 0
    if isempty(XTest) || isempty(yTest)
        % No need to compute this if you're just after the model
        accuracy = [];
        return
    else
        yPredict = predict(Mdl,XTest);
        accuracy = BF_LossFunction(yTest,yPredict,cfnParams.whatLoss,cfnParams.classLabels);
    end
else
    % Test data is mixed through the training data provided using k-fold cross validation
    % Output is the accuracy/loss measure for each fold
    if cfnParams.computePerFold
        % Compute separately for each test fold:
        yPredict = kfoldPredict(Mdl);
        accuracyTestFolds = arrayfun(@(x) BF_LossFunction(yTrain(Mdl.Partition.test(x)),...
                            yPredict(Mdl.Partition.test(x)),cfnParams.whatLoss,...
                            cfnParams.classLabels),1:cfnParams.numFolds);

        % Compute for each training fold:
        accuracyTrainFolds = zeros(cfnParams.numFolds,1);
        for k = 1:cfnParams.numFolds
            XTrain_k = XTrain(Mdl.Partition.training(k),:);
            yTrain_k = yTrain(Mdl.Partition.training(k));
            yPredict_k = predict(Mdl.Trained{k},XTrain_k);
            accuracyTrainFolds(k) = BF_LossFunction(yTrain_k,yPredict_k,cfnParams.whatLoss,...
                            cfnParams.classLabels);
        end

        % Combine train/test accuracies:
        accuracy = {accuracyTrainFolds,accuracyTestFolds};

        % Warning for training fold over-fitting:
        if all(accuracyTrainFolds==100)
            warning('Over-fitting alert (100% in-sample accuracy): consider reducing feature space or regularizing')
        end
    else
        yPredict = kfoldPredict(Mdl);
        % Aggregate all outputs from folds and compute accuracy once from that aggregate:
        accuracy = BF_LossFunction(yTrain,yPredict,cfnParams.whatLoss,cfnParams.classLabels);

        % % Above is equivalent to computing the accuracy in each test fold, then computing the mean
        % % across folds:
        % accuracyTestFolds = arrayfun(@(x) BF_LossFunction(yTrain(Mdl.Partition.test(x)),...
        %                 yPredict(Mdl.Partition.test(x)),cfnParams.whatLoss,...
        %                 cfnParams.classLabels),1:cfnParams.numFolds);
        % accuracy = mean(accuracyTestFolds);
    end
end

end
