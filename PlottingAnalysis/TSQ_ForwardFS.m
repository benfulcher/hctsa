% --------------------------------------------------------------------------
% TSQ_ForwardFS
% --------------------------------------------------------------------------
% 
% Performs greedy forward feature selection for a given classification of the
% data. After selecting the features (using specified training indices), then
% applies the learned classification rule to the training and test sets to get
% training and test classification errors.
% 
% Typical usage uses 'linear' for the criterion (linear classification rates).
% 
%---OUTPUTS:
% ifeat: indices of features selected.
% TestStat: test statistics for all operations.
% TrainErr: training errors for selected features.
% TestErr: test errors for selected features.
% TestClass: classificaiton of the test data.
%
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function [ifeat, TestStat, TrainErr, TestErr, TestClass] = TSQ_ForwardFS(WhatData, ...
                        iTrain,criterion,CrossVal,NumFeatSelect,plotoutput,howzero)

% --------------------------------------------------------------------------
%% Check inputs:
% --------------------------------------------------------------------------
if nargin < 1 || isempty(WhatData)
    error('You must provide data or specify a data source!');
end

if nargin < 2 || isempty(iTrain)
    error('You must specify training indices');
end

if nargin < 3 || isempty(criterion)
    criterion = 'linclass';
    fprintf(1,'Default: Using in-sample linear classification rate\n');
end

if nargin < 4 || isempty(CrossVal)
    CrossVal = 'none';
    fprintf(1,'Default: No cross-validation inside the training data\n');
end

if nargin < 5 || isempty(NumFeatSelect)
    NumFeatSelect = 2; % Stop after two features are selected
end

if nargin < 6 || isempty(plotoutput)
    plotoutput = 1;
end

if nargin < 7 || isempty(howzero)
    % How to deal with operations not improving or hitting zero training
    % error...?
    howzero = 'rand'; % 'rand' (chooses ops at random) ,'NaN' (makes NaNs)
end

% --------------------------------------------------------------------------
%% Load the data
% --------------------------------------------------------------------------
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

% --------------------------------------------------------------------------
%% Set testing indices
% --------------------------------------------------------------------------
% Check whether training indices are logicals
if islogical(iTrain) && length(iTrain)==length(TimeSeries)
    iTrain = find(iTrain); % Convert to indices
end
% Train on iTrain, test on the rest
iTest = setxor((1:length(TimeSeries)),iTrain);
fprintf(1,['We have %u / %u (= %3.1f%%) for training and %u / %u (= %3.1f%%) for ' ...
            'testing.\n'],length(iTrain), length(TimeSeries), ...
            length(iTrain)/length(TimeSeries)*100,length(iTest), ...
                length(TimeSeries),length(iTest)/length(TimeSeries)*100);

% --------------------------------------------------------------------------
%% Set up the classification function
% --------------------------------------------------------------------------
% We get (*) Classify_fn (takes test labels as input): gives the mean classification rate across partitions
%        (*) Classify_fn_label (doesn't take test labels as input): gives labels assigned to test set.

switch criterion
case {'linear','linclass'}
    fprintf(1,'A linear classifier\n');
    Classify_fn_label = @(XTrain,yTrain,Xtest)(classify(Xtest,XTrain,yTrain,'linear'));
    Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest ~= classify(Xtest,XTrain,yTrain,'linear'))/length(ytest);
case 'diaglinear' % Naive Bayes
    fprintf(1,'A Naive Bayes classifier\n');
    Classify_fn_label = @(XTrain,yTrain,Xtest)(classify(Xtest,XTrain,yTrain,'diaglinear'));
    Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest ~= classify(Xtest,XTrain,yTrain,'diaglinear'))/length(ytest);
case {'svm','svmlinear'}
    fprintf(1,'A linear support vector machine\n');
    Classify_fn_label = @(XTrain,yTrain,Xtest) ...
                    svmclassify(svmtrain(XTrain,yTrain, ...
                                'Kernel_Function','linear'),Xtest);
    Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest ~= svmclassify(svmtrain(XTrain,yTrain, ...
                                'Kernel_Function','linear'),Xtest))/length(ytest);
otherwise
    error('Unknown classification method ''%s''',criterion)
end


% --------------------------------------------------------------------------
%% Set up the data partitions for cross-validation
% --------------------------------------------------------------------------

% Define groups as a column vector:
TimeSeriesGroups = [TimeSeries.Group]';

switch CrossVal
case 'kfold'
    fprintf(1,'Using 10-fold stratified cross-validation within the training data\n');
    DataPartitions = cvpartition(TimeSeriesGroups(iTrain),'k',10);
case 'leaveout'
    fprintf(1,'Using leave-one-out cross-validation within the training data\n');
    DataPartitions = cvpartition(TimeSeriesGroups(iTrain),'leaveout');
    fprintf(1,'(Using %u different test sets)\n',DataPartitions.NumTestSets);
case 'none'
    fprintf(1,'No cross-validation performed within the training set\n');
    DataPartitions = [];
otherwise
    error('Unknown cross validation setting ''%s''',CrossVal)
end

% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% Do the feature selection
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------

fprintf(1,['Performing greedy forward feature selection ' ...
                        'using ''%s'' on the training data...\n'],criterion);

% Initialize variables:
ifeat = zeros(NumFeatSelect,1); % Stores indicies of features chosen at each stage
TestStat = zeros(NumFeatSelect,length(Operations)); % Test statistic (classification rate) for all features at each stage
TrainErr = zeros(NumFeatSelect,1); % Training error at each iteration

FS_timer = tic; % start a timer
for j = 1:NumFeatSelect
    % Find the feature that works best in combination with those already chosen:
    for i = 1:length(Operations)
        try
            if isempty(DataPartitions)
                % No cross-validation
                % Compute in-sample training errors and hope we're not overfitting ;-)
                TestStat(j,i) = Classify_fn(TS_DataMat(iTrain,[ifeat(1:j-1);i]),TimeSeriesGroups(iTrain), ...
                                        TS_DataMat(iTrain,[ifeat(1:j-1);i]),TimeSeriesGroups(iTrain));
            else
                % Take mean over cross-validation data partitions specified above:
                TestStat(j,i) = mean(crossval(Classify_fn,TS_DataMat(iTrain,[ifeat(1:j-1);i]), ...
                                    TimeSeriesGroups(iTrain),'partition',DataPartitions));
            end
        catch emsg
            TestStat(j,i) = NaN;
            % keyboard
            fprintf(1,'Error at iteration %u (with feature %u): %s\n',j,i,emsg.identifier)
        end
    end
    
    % --------------------------------------------------------------------------
    % Add the best feature to the list
    % --------------------------------------------------------------------------            
    if all(isnan(TestStat(j,:)))
        ifeat(j) = NaN;
        TrainErr(j) = NaN;
        error('Error selecting feature at iteration %u\n',j)
    else
        TopOps = find(TestStat(j,:)==min(TestStat(j,:)));
        if length(TopOps) > 1 % More than one 'equal best': pick from best at random
            fprintf(1,['Selecting a feature at random from the %u' ...
                ' operations with %4.1f%% error\n'],length(TopOps),min(TestStat(j,:))*100)
            rp = randperm(length(TopOps));
            ifeat(j) = TopOps(rp(1));
        else
            ifeat(j) = TopOps; % Only one best
        end
        TrainErr(j) = min(TestStat(j,:))*100;
        fprintf(1,'Feature %u: %s (%4.2f%% training error)\n',j,Operations(ifeat(j)).Name,TrainErr(j))
    end
    
    % --------------------------------------------------------------------------
    % Evaluate the termination criteria
    % --------------------------------------------------------------------------
    if strcmp(howzero,'NaN') && (j < NumFeatSelect) && (j > 1) && (min(TestStat(j,:))==0)
        % Stop because already at zero error
        % Set ifeat for more features to NaN
        fprintf(1,'Already at perfect classification -- stopping here.\n')
        ifeat(j+1:NumFeatSelect) = NaN;
        for k = j+1:NumFeatSelect
            TestStat(k,:) = NaN; %ones(length(Operations),1)*NaN;
        end
        break
    elseif strcmp(howzero,'NaN') && (j > 1) && min(TestStat(j,:))==min(TestStat(j-1,:))
        % Stop because no improvement from adding this feature
        % Set ifeat for this many features to NaN
        fprintf(1,'No improvement!! -- stopping one feature before here.\n')
        ifeat(j:NumFeatSelect) = NaN; % Remove this one
        for k = j:NumFeatSelect
            % Forget about them:
            TestStat(k,:) = NaN; % ones(length(Operations),1)*NaN; 
        end
        break
    end
end

% Finished selecting features!
fprintf(1,'Feature selection to %u features completed in %s.\n',NumFeatSelect,BF_thetime(toc(FS_timer)))
clear FS_timer;


% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%% Classify the test set using the selected features
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
fprintf(1,'Classifying the test data using the %u selected features...\n',NumFeatSelect);


% Get classification rule from training sets and then apply to test sets
TestErr = zeros(NumFeatSelect,1);
TestCfn = zeros(NumFeatSelect,length(iTest));
for i = 1:NumFeatSelect
    % Get the classification rule from the training data
    try
        % Evaluate the classifier trained on all the training data on the test data:
        TestClass = Classify_fn_label(TS_DataMat(iTrain,ifeat(1:i)),TimeSeriesGroups(iTrain), ...
                                TS_DataMat(iTest,ifeat(1:i)));
        TestErr(i) = mean(TestClass~=TimeSeriesGroups(iTest))*100;
        TestCfn(i,:) = TestClass; % The test classification
    catch
        fprintf(1,'Error classifying train and test data at %u\n',i);
        TrainErr(i) = NaN;
        TestErr(i) = NaN;
        TestClass(i,:) = NaN;
    end
    
    fprintf(1,'Train/Test error rates (%u features): %3.1f%% / %3.1f%%\n', ...
                                i,round(TrainErr(i)),round(TestErr(i)))
end

% % print the top 25
% for i = 1:25
%     disp([Operations(ifeat(i)).Name ' -- ' Operaitons(ifeat(i)).Keywords ' :: ' num2str(TestStat(i))]);
% end

% if plotoutput
%     % plot distributions:
%     figure('color','w'); box('on'); hold on;
%     [f,xi] = ksdensity(TestStat);
%     [f2,xi2] = ksdensity(TestStat2);
%     plot(xi,f,'b');
%     plot(xi2,f2,'r');
%     xlabel('error')
%     legend('individual features','other features with feature 1')
% 
%     % 2d scatter plot
%     TSQ_plot_2d(kwgs,gi,ifeat,norcl);
% end


end