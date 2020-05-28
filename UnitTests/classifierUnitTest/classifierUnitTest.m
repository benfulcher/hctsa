% Test the classifier yields the correct accuracy

clear;

rng(1);

%% Configure

numFeatures = 3;

prefix = 'uTest_';
datfile = [prefix 'ts.mat'];
groups = {'medical','sound'};

classifierFilename = [prefix 'classifier.mat'];

%% Setup

hctsa_dir = '../../';
if ~exist('TS_Compute','file')
  cwd = pwd;
  cd(hctsa_dir);
  startup
  cd(cwd);
end

whichOps = randsample(7700,numFeatures); % Choose 3 random operations to use for classification/prediction

TS_Init(datfile,[],[],0,[prefix 'HCTSA.mat']);

TS_Compute(1,[],whichOps,[],[prefix 'HCTSA.mat'],0);

%% For classification

matfile = [prefix 'HCTSA.mat'];

TS_LabelGroups(matfile,groups,true,true);

if exist([prefix 'HCTSA_filtered.mat'],'file')
    matfile = [prefix 'HCTSA_filtered.mat'];
end

matfileNormed = [matfile(1:end-4) '_N.mat'];

% Filter the data to remove features that are null
TS_Normalize('none',[0,1],matfile,true);

% Load normalized data in a structure:
unnormalizedData = load(matfileNormed);

% Set how to normalize the data:
whatNormalization = 'zscore'; % 'zscore', 'scaledRobustSigmoid'

% Normalize the data (for joint classifier):
TS_Normalize(whatNormalization,[0,1],matfile,true);

% Load normalized data in a structure:
normalizedData = load(matfileNormed);

%-------------------------------------------------------------------------------
%% Run classifiers:

whatClassifier = 'svm_linear';
TS_Classify(normalizedData,whatClassifier,'numPCs',0,'numNulls',0,...
            'classifierFilename',classifierFilename);

numFeatures = 40; % number of features to include in the pairwise correlation plot
numFeaturesDistr = 32; % number of features to show class distributions for
whatStatistic = 'fast_linear'; % classification statistic

TS_TopFeatures(unnormalizedData,whatStatistic,'numFeatures',numFeatures,...
                'numFeaturesDistr',numFeaturesDistr,...
                'whatPlots',{'histogram','distributions','cluster'},...
                'classifierFilename',classifierFilename);
            
%% Run predictions (should match classifier accy):
% Note we can use just acc (not balancedAcc) since the classes are balanced

topFeatureTab = TS_Predict(unnormalizedData.TimeSeries.Data,...
                    unnormalizedData.TimeSeries.Name,...
                    classifierFilename,...
                    'classifierType','topFeature');

jointTab = TS_Predict(normalizedData.TimeSeries.Data,...
                    normalizedData.TimeSeries.Name,...
                    classifierFilename,...
                    'classifierType','allFeatures');
                
load(classifierFilename);

topFeatureAcc = 100*mean(topFeatureTab.predictGroups == normalizedData.TimeSeries.Group);
jointAcc = 100*mean(jointTab.predictGroups == normalizedData.TimeSeries.Group);
expTopFeatureAcc = featureClassifier.Accuracy;
expJointAcc = jointClassifier.Accuracy;

%% Check the output
fprintf('\n++ Unit Tests ++\n');
fprintf('Checking joint-feature classifier prediction loss matches training loss...')
if jointAcc == expJointAcc
    fprintf(' passed.\n');
else
    fprintf(' failed (expected: %.3f, received: %.3f).\n', expJointAcc, jointAcc);
end
fprintf('Checking top-feature classifier prediction loss matches training loss...')
if topFeatureAcc == expTopFeatureAcc
    fprintf(' passed.\n');
else
    fprintf(' failed (expected: %.3f, received: %.3f).\n', expTopFeatureAcc, topFeatureAcc);
end

%% clean up

out = input(sprintf('Remove periphery files (%s:%s:%s)? [yn] ',matfile,matfileNormed,classifierFilename),'s');
if out == 'y'
    delete(matfile);
    delete(matfileNormed);
    delete(classifierFilename);
end