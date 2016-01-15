function [ifeat, testStat, trainErr, testErr, TestClass] = TS_ForwardFS(whatData,iTrain,criterion,crossVal,numFeatSelect,howzero)
% TS_ForwardFS  Greedy forward feature selection
%
% Uses the sequentialfs function from Matlab's Statistics Toolbox.
% After selecting the features (using specified training indices), then
% applies the learned classification rule to the training and test sets to get
% training and test classification errors.
%
% Typical usage uses 'linear' for the criterion (linear classification rates),
% or 'diaglinear' when there are highly correlated features in the data matrix.
%
% NOTE: This function requires a training portion to be specified. If there
%       is no obvious training portion of your data, this should be specified
%       at random (e.g., stratified to represent groups in your data) and
%       repeated to get an estimate of the out-of-sample classification that is
%       not dependent on the particular data partition.
%
%---INPUTS:
% whatData: the data to load in (cf. TS_LoadData)
% iTrain: training indices
% criterion: what criterion on which to evaluate the quality of a feature set
%
%---OUTPUTS:
% ifeat: indices of features selected.
% testStat: test statistics for all operations.
% trainErr: training errors for selected features.
% testErr: test errors for selected features.
% TestClass: classificaiton of the test data.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
end

if nargin < 2
    iTrain = [];
end

if nargin < 3 || isempty(criterion)
    criterion = 'linclass';
    fprintf(1,'Default: Using in-sample linear classification rate\n');
end

if nargin < 4 || isempty(crossVal)
    crossVal = 'none';
    fprintf(1,'Default: No cross-validation inside the training data\n');
end

if nargin < 5 || isempty(numFeatSelect)
    numFeatSelect = 2; % Stop after two features are selected
end

if nargin < 6 || isempty(howzero)
    % How to deal with operations not improving or hitting zero training
    % error...?
    % (*) 'rand' (chooses ops at random)
    % (*) 'NaN' (makes NaNs)
    howzero = 'rand';
end

% ------------------------------------------------------------------------------
%% Load the data
% ------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);

% Retrieve group names also:
fileVarsStruct = whos('-file',whatDataFile);
fileVars = {fileVarsStruct.name};
if ismember('groupNames',fileVars)
    groupNames = load(whatDataFile,'groupNames');
    groupNames = groupNames.groupNames;
else
    groupNames = {};
end

%-------------------------------------------------------------------------------
% Set training indices (if empty)
%-------------------------------------------------------------------------------
if isempty(iTrain)
    warning('Using a random stratified 80% sample of the data for training...')
    cvp = cvpartition([TimeSeries.Group],'holdout',0.2);
    iTrain = cvp.training;
end

% ------------------------------------------------------------------------------
%% Set testing indices
% ------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
%% Set up the classification function
% ------------------------------------------------------------------------------
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


% ------------------------------------------------------------------------------
%% Set up the data partitions for cross-validation
% ------------------------------------------------------------------------------

% Define groups as a column vector:
TimeSeriesGroups = [TimeSeries.Group]';

switch crossVal
case 'kfold'
    fprintf(1,'Using 10-fold stratified cross-validation within the training data\n');
    dataPartitions = cvpartition(TimeSeriesGroups(iTrain),'k',10);
case 'leaveout'
    fprintf(1,'Using leave-one-out cross-validation within the training data\n');
    dataPartitions = cvpartition(TimeSeriesGroups(iTrain),'leaveout');
    fprintf(1,'(Using %u different test sets)\n',dataPartitions.NumTestSets);
case 'none'
    fprintf(1,'No cross-validation performed within the training set\n');
    dataPartitions = [];
otherwise
    error('Unknown cross validation setting ''%s''',crossVal)
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
%% Do the feature selection
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

fprintf(1,['Performing greedy forward feature selection ' ...
                        'using ''%s'' on the training data...\n'],criterion);

% Initialize variables:
ifeat = zeros(numFeatSelect,1); % Stores indicies of features chosen at each stage
testStat = zeros(numFeatSelect,length(Operations)); % Test statistic (classification rate) for all features at each stage
trainErr = zeros(numFeatSelect,1); % Training error at each iteration

FS_timer = tic; % start a timer
for j = 1:numFeatSelect
    % Find the feature that works best in combination with those already chosen:
    for i = 1:length(Operations)
        try
            if isempty(dataPartitions)
                % No cross-validation
                % Compute in-sample training errors and hope we're not overfitting ;-)
                testStat(j,i) = Classify_fn(TS_DataMat(iTrain,[ifeat(1:j-1);i]),TimeSeriesGroups(iTrain), ...
                                        TS_DataMat(iTrain,[ifeat(1:j-1);i]),TimeSeriesGroups(iTrain));
            else
                % Take mean over cross-validation data partitions specified above:
                testStat(j,i) = mean(crossval(Classify_fn,TS_DataMat(iTrain,[ifeat(1:j-1);i]), ...
                                    TimeSeriesGroups(iTrain),'partition',dataPartitions));
            end
        catch emsg
            testStat(j,i) = NaN;
            fprintf(1,'Error at iteration %u (with feature %u): %s\n',j,i,emsg.identifier)
        end
    end

    % ------------------------------------------------------------------------------
    % Add the best feature to the list
    % ------------------------------------------------------------------------------
    if all(isnan(testStat(j,:)))
        ifeat(j) = NaN;
        trainErr(j) = NaN;
        error('Error selecting feature at iteration %u\n',j)
    else
        TopOps = find(testStat(j,:)==min(testStat(j,:)));
        if length(TopOps) > 1 % More than one 'equal best': pick from best at random
            fprintf(1,['Selecting a feature at random from the %u' ...
                ' operations with %4.1f%% error\n'],length(TopOps),min(testStat(j,:))*100)
            rp = randperm(length(TopOps));
            ifeat(j) = TopOps(rp(1));
        else
            ifeat(j) = TopOps; % Only one best
        end
        trainErr(j) = min(testStat(j,:))*100;
        fprintf(1,'Feature %u: %s (%4.2f%% training error)\n',j,Operations(ifeat(j)).Name,trainErr(j))
    end

    % ------------------------------------------------------------------------------
    % Evaluate the termination criteria
    % ------------------------------------------------------------------------------
    if strcmp(howzero,'NaN') && (j < numFeatSelect) && (j > 1) && (min(testStat(j,:))==0)
        % Stop because already at zero error
        % Set ifeat for more features to NaN
        fprintf(1,'Already at perfect classification -- stopping here.\n')
        ifeat(j+1:numFeatSelect) = NaN;
        for k = j+1:numFeatSelect
            testStat(k,:) = NaN; %ones(length(Operations),1)*NaN;
        end
        break
    elseif strcmp(howzero,'NaN') && (j > 1) && min(testStat(j,:))==min(testStat(j-1,:))
        % Stop because no improvement from adding this feature
        % Set ifeat for this many features to NaN
        fprintf(1,'No improvement!! -- stopping one feature before here.\n')
        ifeat(j:numFeatSelect) = NaN; % Remove this one
        for k = j:numFeatSelect
            % Forget about them:
            testStat(k,:) = NaN; % ones(length(Operations),1)*NaN;
        end
        break
    end
end

% Finished selecting features!
fprintf(1,'Feature selection to %u features completed in %s.\n',...
            numFeatSelect,BF_thetime(toc(FS_timer)))
clear FS_timer;

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
%% Classify the test set using the selected features
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
fprintf(1,'Classifying the test data using the %u selected features...\n',numFeatSelect);


% Get classification rule from training sets and then apply to test sets
testErr = zeros(numFeatSelect,1);
testCfn = zeros(numFeatSelect,length(iTest));
for i = 1:numFeatSelect
    % Get the classification rule from the training data
    try
        % Evaluate the classifier trained on all the training data on the test data:
        TestClass = Classify_fn_label(TS_DataMat(iTrain,ifeat(1:i)),TimeSeriesGroups(iTrain), ...
                                TS_DataMat(iTest,ifeat(1:i)));
        testErr(i) = mean(TestClass~=TimeSeriesGroups(iTest))*100;
        testCfn(i,:) = TestClass; % The test classification
    catch
        fprintf(1,'Error classifying train and test data at %u\n',i);
        trainErr(i) = NaN;
        testErr(i) = NaN;
        TestClass(i,:) = NaN;
    end

    fprintf(1,'Train/Test error rates (%u features): %3.1f%% / %3.1f%%\n', ...
                                i,round(trainErr(i)),round(testErr(i)))
end

%-------------------------------------------------------------------------------
% Plot in a 2-D space
%-------------------------------------------------------------------------------
% Set feature labels:
featureLabels = cell(2,1);
for i = 1:2
    featureLabels{i} = sprintf('%s (%.2f%%/%.2f%%)',Operations(ifeat(i)).Name,trainErr(i),testErr(i));
end

TS_plot_2d(TS_DataMat(:,ifeat(1:2)),TimeSeries,featureLabels,groupNames)


end
