function TS_classify(whatData,whatClassifier,doPCs,seedReset)
% TS_classify   Classify groups in the data using all features (and PCs)
%
% This function uses a classifier to learn group labels assigned to time series
% in the dataset.
%
%---USAGE:
% TS_classify;
%
%---INPUTS:
% whatData, the hctsa data to use (input to TS_LoadData)
% whatClassifier, the classification method to use ('svm', 'knn', 'linear')
% doPCs, (binary) whether to investigate the behavior of different numbers of PCs of the
%        data matrix (default: 1)
% seedReset, input to BF_ResetSeed, specifying whether (and how) to reset the
%               random seed (for reproducible results from cross-validation)
%
%---OUTPUTS:
% Text output on classification rate using all features, and if doPCs = 1, also
% shows how this varies as a function of reduced PCs (text and as an output plot)

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
% Check Inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
if nargin < 2
    whatClassifier = 'svm';
    % 'svm', 'discriminant', 'knn'
end
if nargin < 3
    doPCs = 1;
end
if nargin < 4
    seedReset = 'default';
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

% Check that group labels have been assigned
if ~isfield(TimeSeries,'Group')
    error('Group labels not assigned to time series. Use TS_LabelGroups.');
end
timeSeriesGroup = [TimeSeries.Group]; % Use group form
numClasses = length(unique(timeSeriesGroup));
numFeatures = length(Operations);
numTimeSeries = length(TimeSeries);

%-------------------------------------------------------------------------------
% Set up the classification model
%-------------------------------------------------------------------------------
switch whatClassifier
case 'svm'
    % Linear SVM:
    cfnModel = templateSVM('Standardize',1,'KernelFunction','linear');
case 'knn'
    % k-NN (k=3) classification:
    cfnModel = templateKNN('NumNeighbors',3,'Distance','euclidean');
case {'discriminant','linear'}
    % Linear discriminant analysis:
    cfnModel = templateDiscriminant('DiscrimType','linear');
    % could also be 'naivebayes', 'tree', ensemble methods
otherwise
    error('Unknown classification model, ''%s''',whatClassifier);
end

%-------------------------------------------------------------------------------
% Fit the model using k-fold cross validation:
%-------------------------------------------------------------------------------

% Reset the seed?
BF_ResetSeed(seedReset);

% Set the number of folds for k-fold cross validation using a heuristic
% (for small datasets with fewer than 10 examples per class):
pointsPerClass = numTimeSeries/numClasses;
if pointsPerClass < 5
    numFolds = 2;
elseif pointsPerClass < 10
    numFolds = 5;
else
    numFolds = 10;
end

% Fit the classification model to the dataset (for each cross-validation fold)
CVcfnModel = fitcecoc(TS_DataMat,timeSeriesGroup,'Learners',cfnModel,'KFold',numFolds);

% Get the misclassification rate from each fold:
foldLosses = 100*(1 - kfoldLoss(CVcfnModel,'Mode','individual'));

fprintf(1,['\nClassification rate (%u-class) using %u-fold %s classification with %u' ...
                 ' features:\n%.3f +/- %.3f%%\n\n'],...
                            numClasses,...
                            numFolds,...
                            whatClassifier,...
                            numFeatures,...
                            mean(foldLosses),...
                            std(foldLosses))


% f = figure('color','w');
% histogram(foldLosses*100)
% xlim([0,100]);

%-------------------------------------------------------------------------------
% Compare performance of reduced PCs from the data matrix:
%-------------------------------------------------------------------------------
if doPCs
    numPCs = 10;

    % Compute top 10 PCs of the data matrix:
    fprintf('Computing top %u PCs...',numPCs)
    [pcCoeff, pcScore, latent, ~, percVar] = pca(zscore(TS_DataMat),'NumComponents',numPCs);
    fprintf(' Done.\n')
    numPCs = min(10,size(pcScore,2)); % sometimes lower than attempted 10

    % Compute cumulative performance of PCs:
    PC_cfn = cell(numPCs,1);
    cfnRate = zeros(numPCs,2);
    fprintf('Computing classification rates keeping top 1--%u PCs...\n',numPCs)
    for i = 1:numPCs
        PC_cfn{i} = fitcecoc(pcScore(:,1:i),timeSeriesGroup,'Learners',cfnModel,'KFold',numFolds);
        losses = 1-kfoldLoss(PC_cfn{i},'Mode','individual');
        cfnRate(i,1) = mean(losses)*100;
        cfnRate(i,2) = std(losses)*100;
        fprintf(1,'%u PCs:   %.3f +/- %.3f%%\n',i,cfnRate(i,1),cfnRate(i,2))
    end

    plotColors = BF_getcmap('spectral',3,1);

    f = figure('color','w'); hold on
    plot([1,numPCs],ones(2,1)*mean(foldLosses),'--','color',plotColors{3})
    plot(1:numPCs,cfnRate(:,1),'o-k')
    legend(sprintf('All %u features (%.1f%%)',numFeatures,mean(foldLosses)),...
                sprintf('PCs (%.1f--%.1f%%)',min(cfnRate(:,1)),max(cfnRate(:,1))))
    plot(1:numPCs,cfnRate(:,1)+cfnRate(:,2),':k')
    plot(1:numPCs,cfnRate(:,1)-cfnRate(:,2),':k')

    xlabel('Number of PCs');
    ylabel('Classification accuracy (%)')

    titleText = sprintf('Classification rate (%u-class) using %u-fold %s classification',...
                                numClasses,...
                                numFolds,...
                                whatClassifier);
    title(titleText)

end

end
