% TSQ_ForwardFS
% 
% Performs greedy forward feature selection for a given classification of the
% data. After selecting the features (using specified training indices), then
% applies the learned classification rule to the training and test sets to get
% training and test classification errors.
% 
% Typical usage uses 'linclass' or 'linclasscv' for the criterion.
% 
%--OUTPUTS:
%-ifeat: indices of features selected
%-teststatout: test statistics for all operations.
%
%-------HISTORY----------
% uses BioInformatics toolbox functions
% Finds the best combination of 2 features for the classification task
% Ben Fulcher 23/9/2010
% Ben Fulcher 14/10/2010: added subset
% Ben Fulcher 20/10/2010: added nbest (return this many best features)
% 
% Some options:
% 'ttest' (default) -- Absolute value two-sample t-test with pooled variance
%                      estimate.
% 'entropy' -- Relative entropy, also known as Kullback-Leibler distance or
%              divergence.
% 'bhattacharyya' -- Minimum attainable classification error or Chernoff bound.
% 'roc' -- Area between the empirical receiver operating characteristic (ROC)
%          curve and the random classifier slope.
% 'wilcoxon' -- Absolute value of the u-statistic of a two-sample unpaired
%               Wilcoxon test, also known as Mann-Whitney.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function [ifeat, teststat, TrainErr, TestErr, TestClass] = TSQ_ForwardFS(WhatData, ...
                            iTrain,criterion,randomize,plotoutput,NumFeatSelect,howzero,knn)

%% Inputs
if nargin < 1 || isempty(WhatData)
    error('You must provide data or specify a data source!');
end

if nargin < 2 || isempty(iTrain)
    error('You must specify training indices');
end

if nargin < 3 || isempty(criterion)
    criterion = 'roc';
    fprintf(1,'Default: Using roc criterion\n');
end

if nargin < 4 || isempty(randomize)
    randomize = 0;
end

if nargin < 5 || isempty(plotoutput)
    plotoutput = 1;
end

if nargin < 6 || isempty(NumFeatSelect)
    NumFeatSelect = 2; % Stop after two features are selected
end

if nargin < 7 || isempty(howzero)
    % How to deal with operations not improving or hitting zero training
    % error...?
    howzero = 'rand'; % 'rand' (chooses ops at random) ,'NaN' (makes NaNs)
end

if nargin < 8
    knn = 1; % knn(1) by default
end

%% Load the data
% If specified a string: 'norm' or 'cl' will retrieve
% Otherwise must be a structure with fields 'TimeSeries', 'Operations' and
% 'TS_DataMat'
if ischar(WhatData)
    switch WhatData
    case 'norm'
        TheDataFile = 'HCTSA_N.mat';
    case 'cl'
        TheDataFile = 'HCTSA_cl.mat';
    end
    fprintf(1,'Loading data from %s...',TheDataFile);
    load(TheDataFile,'TS_DataMat','Operations','TimeSeries');
    fprintf(1,' Loaded.\n');
elseif isstruct(WhatData)
    % Already loaded and given here as a structure
    TS_DataMat = WhatData.TS_DataMat;
    Operations = WhatData.Operations;
    TimeSeries = WhatData.TimeSeries;
    fprintf(1,'Data provided, adapted successfully.\n');
end

% Set testing indices
% Train on iTrain, test on the rest
iTest = setxor((1:length(TimeSeries)),iTrain);
fprintf(1,'We have %u / %u training data items (= %4.1f%%)\n',length(iTrain), ...
                    length(TimeSeries),length(iTrain)/length(TimeSeries)*100);
fprintf(1,'We have %u / %u testing data items (= %4.1f%%)\n',length(iTest), ...
                    length(TimeSeries),length(iTest)/length(TimeSeries)*100);

% Use random numbers instead of the actual values in the data matrix
if randomize
    fprintf(1,'Using random numbers instead of the actual information in the data matrix\n');
    TS_DataMat = randn(size(TS_DataMat));
end

%% Run the algorithm
TimeSeriesGroups = [TimeSeries.Group]'; % one of my functions to convert gi to group form

FS_timer = tic; % start a timer

fprintf(1,'Performing feature selection using ''%s''...\n',criterion);

switch criterion
    case {'linclasscv','linclass'}
% if strcmp(criterion,'linclasscv') || strcmp(criterion,'linclass') % quantify linear classification for each
% F_linclass = @(XT,yT,Xt,yt) sum(yt~=classify(Xt,XT,yT,'linear'))/length(yt);
        classf = @(XTRAIN,ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain,'linear'));
        if strcmp(criterion,'linclasscv')
            % 10-fold stratified cross-validation WITHIN THE TRAINING DATA
            DataPartitions = cvpartition(TimeSeriesGroups(iTrain),'k',10);
        end

        ifeat = zeros(NumFeatSelect,1); % stores indicies of features chosen at each stage
        teststat = cell(NumFeatSelect,1); % cross-validation classification rates for all features at each stage
        for j = 1:NumFeatSelect
            % Find the feature that works best in combination with those already chosen:
            teststat{j} = zeros(length(Operations),1);
            for i = 1:length(Operations)
                try
                    if strcmp(criterion,'linclass')
                        % in-sample errors
                        [~,err] = classify(TS_DataMat(iTrain,[ifeat(1:j-1);i]),TS_DataMat(iTrain,[ifeat(1:j-1);i]), ...
                                                TimeSeriesGroups(iTrain),'linear');
                        teststat{j}(i) = err;
                    else
            %             mcrs = crossval(TS_DataMat_linclass,TS_DataMat(:,i),TimeSeriesGroups,'partition',DataPartitions);
                        teststat{j}(i) = crossval('mcr',TS_DataMat(iTrain,[ifeat(1:j-1);i]),TimeSeriesGroups(iTrain), ...
                                                    'predfun',classf,'partition',DataPartitions);
                        % This code with 'mcr' is the same as mean(mcrs)
                    end
                catch emsg
                    teststat{j}(i) = NaN;
    %                 disp(emsg.identifier)
        %             figure('color','w');hold on
        %             plot_ks(TS_DataMat(gi{1},i),[1,0,0],0) % plot the first group
        %             plot_ks(TS_DataMat(gi{2},i),[0,0,1],0) % plot the second group
        %             title(Operations(i).Name);
        %             keyboard
                end
            end
            if all(isnan(teststat{j}))
                ifeat(j) = NaN;
                fprintf(1,'Error selecting feature at iteration %u\n',j)
                keyboard
            else
                tops = find(teststat{j}==min(teststat{j}));
                if length(tops) > 1 % more than one 'equal best': pick from best at random
                    fprintf(1,['Selecting a feature at random from the %u' ...
                        ' operations with %4.1f%% error\n'],length(tops),min(teststat{j})*100)
                    rp = randperm(length(tops));
                    ifeat(j) = tops(rp(1));
                else
                    ifeat(j) = tops; % only one best
                end
                fprintf(1,'Feature %u: %s (%4.2f%%)\n',j,Operations(ifeat(j)).Name,min(teststat{j})*100)
            end
            
            if strcmp(howzero,'NaN') && (j < NumFeatSelect) && (j > 1) && (min(teststat{j})==0)
                % stop because already at zero error
                % set ifeat for more features to NaN
                fprintf(1,'Already at perfect classification -- stopping here.\n')
                ifeat(j+1:NumFeatSelect) = NaN;
                for k = j+1:NumFeatSelect
                    teststat{k} = ones(length(Operations),1)*NaN;
                end
                break
            elseif strcmp(howzero,'NaN') && (j > 1) && min(teststat{j})==min(teststat{j-1})
                % Stop because no improvement from adding this feature
                % Set ifeat for this many features to NaN
                fprintf(1,'No improvement!! -- stopping one feature before here.\n')
                ifeat(j:NumFeatSelect) = NaN; % remove this one
                for k = j:NumFeatSelect
                    teststat{k} = ones(length(Operations),1)*NaN; % Forget about them
                end
                break
            end
        end
        
%     case 'svm_matlab'
%         disp('BioInf Matlab toolbox svm')
%         
%         ifeat = zeros(NumFeatSelect,1); % stores indicies of features chosen at each stage
%         teststat = cell(NumFeatSelect,1); % cross-validation classification rates for all features at each stage
%         for j = 1:NumFeatSelect
%             teststat{j} = zeros(length(Operations),1);
%             for i = 1:length(Operations)
%                 % train linear SVM model
%                 svmStruct = svmtrain(TS_DataMat(:,[ifeat(1:j-1);i]),TimeSeriesGroups,'Kernel_Function','linear');
%                 % classify the same data with the SVM trained on the data
%                 predictedclasses = svmclassify(svmStruct,TS_DataMat(:,[ifeat(1:j-1);i]));
%                 teststat{j}(i) = 1 - sum(predictedclasses==TimeSeriesGroups)/length(TimeSeriesGroups);
%             end
%             ifeat(j) = find(teststat{j}==min(teststat{j}),1);
%             fprintf(1,'Feature %u = %s (%4.1f%%)\n',j,Operations(ifeat(j)).Name,min(teststat{j})*100)
%         end
%         SVMStruct = svmtrain(TS_DataMatsub(itrain,:),TimeSeriesGroups(itrain),'Kernel_Function','linear');
%         
%     case {'knn_matlab','knn'}
% %         disp('BEN KNN(1). This will always give zero error cos it''s knn(1) in-sample!')
% %         disp('KNN(1)')
% %         k = 1;
%         
%         disp(['KNN(' num2str(knn) ')'])
%         clear opts; opts.distance = 'Euclidean'; opts.traintrain = 1;
%         ifeat = zeros(NumFeatSelect,1); % stores indicies of features chosen at each stage
%         teststat = cell(NumFeatSelect,1); % cross-validation classification rates for all features at each stage
%         for j = 1:NumFeatSelect
%             teststat{j} = zeros(size(TS_DataMat,2),1);
%             for i = 1:size(TS_DataMat,2)
%                 [~,err] = benknn(TS_DataMat(:,[ifeat(1:j-1);i]),TS_DataMat(:,[ifeat(1:j-1);i]),knn,TimeSeriesGroups,TimeSeriesGroups,opts); % classifies in-sample
%                 teststat{j}(i) = err;
%             end
%             tops = find(teststat{j}==min(teststat{j}));
%             if length(tops) > 1 % more than one 'equal best': pick from best at random
%                 disp(['Picking a feature at random from the ' num2str(length(tops)) ...
%                     ' operations with ' num2str(min(teststat{j})*100) '% error'])
%                 rp = randperm(length(tops));
%                 ifeat(j) = tops(rp(1));
%             else
%                 ifeat(j) = tops; % only one best
%             end
%             disp(['feature ' num2str(j) ' = ' Operations(ifeat(j)).Name ' (' num2str(min(teststat{j})*100) ' %)'])
%         end
%         
%     case 'knncv'
%         % uses matlab KNN(3) with 5-fold cross-validation
%         disp(['KNN(' num2str(knn) ')'])
%         kfolds = 5;
%         knn = 1;
%         DataPartitions = cvpartition(TimeSeriesGroups,'kfold',kfolds); % specify statified k-fold crossvalidation
%         ifeat = zeros(NumFeatSelect,1); % stores indicies of features chosen at each stage
%         teststat = cell(NumFeatSelect,1); % cross-validation classification rates for all features at each stage
%         for j = 1:NumFeatSelect
%             teststat{j} = zeros(size(TS_DataMat,2),1);
%             for i = 1:size(TS_DataMat,2)
%                 errs = zeros(kfolds,1);
%                 for k = 1:kfolds
%                     itrain = training(DataPartitions,k); % indicies for training
%                     itest = test(DataPartitions,k); % indicies for testing
%                     [~,errs(k)] = benknn(TS_DataMat(itrain,[ifeat(1:j-1);i]),TS_DataMat(itest,[ifeat(1:j-1);i]),knn,TimeSeriesGroups(itrain),TimeSeriesGroups(itest)); % classifies in-sample
%                 end
%                 teststat{j}(i) = mean(errs); % return mean over the cross-validation folds
%             end
%             ifeat(j) = find(teststat{j}==min(teststat{j}),1);
%             disp(['feature ' num2str(j) ' = ' Operations(ifeat(j)).Name ' (' num2str(min(teststat{j})*100) ' %)'])
%         end
%         
%     case {'knn_spider','svm'}
%         % define the spider model
%         if strcmp(criterion,'knn')
%             a = SPIDER_getmemodel('knn',3); % define a knn(3) model
%         else
%             a = SPIDER_getmemodel('svm',{{'linear'}}); % define a svm (linear) model
%         end
%         
%         % initialize the variables
%         ifeat = zeros(NumFeatSelect,1); % stores indicies of features chosen at each stage
%         teststat = cell(NumFeatSelect,1); % cross-validation classification rates for all features at each stage
%         
%         % calculate losses across all features (in combination with those
%         %                                   already chosen)
%         for j = 1:NumFeatSelect
%             % Find the feature that works best in combination with those already chosen:
%             teststat{j} = zeros(size(TS_DataMat,2),1);
%             for i = 1:size(TS_DataMat,2)        
%                 dtrain = makeitdata(TS_DataMat(:,[ifeat(1:j-1);i]),gi);
%                 try
%                     [~,a] = train(a,dtrain); % train classification model on training data
%                     lossme = loss(test(a,dtrain),'class_loss'); % get in-sample loss
%                     teststat{j}(i) = lossme.Y;
%                 catch
%                     teststat{j}(i) = NaN;
%                 end
%             end
%             ifeat(j) = find(teststat{j}==min(teststat{j}),1);
%             disp(['feature ' num2str(j) ' = ' Operations(ifeat(j)).Name ' (' num2str(min(teststat{j})*100) ' %)'])
%         end
%         
%     otherwise
%         fprintf(1,'Using ''rankfeatures'', which is not a custom classification method\n')
%         [ifeat,teststat] = rankfeatures(TS_DataMat',TimeSeriesGroups,'criterion',criterion);
%         [teststat,ix] = sort(teststat,'descend');
%         ifeat = ifeat(ix);
end

% Finished selecting features!
fprintf(1,'Feature selection to %u features complete! Took %s.\n',NumFeatSelect,BF_thetime(toc(FS_timer)))
clear FS_timer;
% Now we have teststat and ifeat

% Format teststat for output
if iscell(teststat)
    teststat = cellfun(@(x)min(x),teststat);
end


%% Now classify the test set using the selected features
fprintf(1,'Classifying the test data using the %u selected features...\n',NumFeatSelect);
switch criterion
    case {'linclass','linclasscv'}
        % Get classification rule from training sets and then apply to test sets
        trainerr = zeros(length(ifeat),1);
        testerr = zeros(length(ifeat),1);
        testclass = zeros(length(ifeat),length(iTest));
        for i = 1:length(ifeat)
            % Get the classification rule from the training data
            try
                [class,err,~,~,coeff] = classify(TS_DataMat(iTest,ifeat(1:i)),TS_DataMat(iTrain,ifeat(1:i)), ...
                                                TimeSeriesGroups(iTrain),'linear');
            
                TrainErr(i) = err*100;
                TestErr(i) = (1 - sum(class==TimeSeriesGroups(iTest))/length(iTest))*100;
                TestClass(i,:) = class; % The test classification
                %(class==gigte); % was this object correctly classified by this feature?
            catch
                fprintf(1,'Error classifying train and test data at %u\n',i);
                TrainErr(i) = NaN;
                TestErr(i) = NaN;
                TestClass(i,:) = NaN;
            end
% %             threshold(i) = -K/L;
            
            fprintf(1,'Train/Test error (%u features): %4.1f%% / %4.1f%%\n', ...
                                        i,round(TrainErr(i)),round(TestErr(i)))
        end

end

% % print the top 25
% for i = 1:25
%     disp([Operations(ifeat(i)).Name ' -- ' Operaitons(ifeat(i)).Keywords ' :: ' num2str(teststat(i))]);
% end

% if plotoutput
%     % plot distributions:
%     figure('color','w'); box('on'); hold on;
%     [f,xi] = ksdensity(teststat);
%     [f2,xi2] = ksdensity(teststat2);
%     plot(xi,f,'b');
%     plot(xi2,f2,'r');
%     xlabel('error')
%     legend('individual features','other features with feature 1')
% 
%     % 2d scatter plot
%     TSQ_plot_2d(kwgs,gi,ifeat,norcl);
% end


end