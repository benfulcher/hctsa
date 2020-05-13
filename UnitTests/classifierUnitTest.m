% Test the classifier yields the correct accuracy

clear;

rng(1);

%% Configure

prefix = './uTest_';
datfile = [prefix 'ts.mat'];
groups = {'medical','sound'};
whichOps = randsample(7700,3); % Choose 3 random operations to use for classification/prediction

classifierFilename = [prefix 'classifier.mat'];

%% Setup

hctsa_dir = '../';
if ~exist('TS_compute','file')
  cwd = pwd;
  cd(hctsa_dir);
  startup
  cd(cwd);
end

TS_init(datfile,[],[],1,[prefix 'HCTSA.mat']);

TS_compute(1,[],whichOps,[],[prefix 'HCTSA.mat']);

%% For classification

matfile = [prefix 'HCTSA.mat'];

TS_LabelGroups(matfile,groups,true,true);

% Set how to normalize the data:
whatNormalization = 'zscore'; % 'zscore', 'scaledRobustSigmoid'

if exist([prefix 'HCTSA_filtered.mat'],'file')
    matfile = [prefix 'HCTSA_filtered.mat'];
end

% Normalize the data, filtering out features with any special values:
TS_normalize(whatNormalization,[0,1],matfile,true);

% Load normalized data in a structure:
normalizedData = load([matfile(1:end-4) '_N.mat']);

%-------------------------------------------------------------------------------
%% Run classifiers:

whatClassifier = 'svm_linear';
TS_classify(normalizedData,whatClassifier,'numPCs',0,'numNulls',0,...
            'classifierFilename',classifierFilename);

numFeatures = 40; % number of features to include in the pairwise correlation plot
numFeaturesDistr = 32; % number of features to show class distributions for
whatStatistic = 'fast_linear'; % classification statistic

TS_TopFeatures(normalizedData,whatStatistic,'numFeatures',numFeatures,...
                'numFeaturesDistr',numFeaturesDistr,...
                'whatPlots',{'histogram','distributions','cluster'},...
                'classifierFilename',classifierFilename);
            
%% Run predictions (should match classifier accy):
% Note we can use just acc (not balancedAcc) since the classes are balanced

topFeatureTab = TS_predict(normalizedData.TimeSeries.Data,...
                    normalizedData.TimeSeries.Name,...
                    classifierFilename,...
                    'classifierType','topFeature',...
                    'predictionFilename','uTest_pred_top.mat');

jointTab = TS_predict(normalizedData.TimeSeries.Data,...
                    normalizedData.TimeSeries.Name,...
                    classifierFilename,...
                    'classifierType','allFeatures',...
                    'predictionFilename','uTest_pred_joint.mat');
                
load(classifierFilename);

topFeatureAcc = 100*mean(topFeatureTab.predictGroups == normalizedData.TimeSeries.Group);
jointAcc = 100*mean(jointTab.predictGroups == normalizedData.TimeSeries.Group);
expTopFeatureAcc = featureClassifier.Accuracy;
expJointAcc = jointClassifier.Accuracy;

fprintf('==============================================================\n');
fprintf('JointFeature classifier acc: %.3f (expected: %.3f)\n',jointAcc,expJointAcc);
fprintf('TopFeature classifier acc: %.3f (expected: %.3f)\n',topFeatureAcc,expTopFeatureAcc);