function outputFileName = TS_Normalize(normFunction,filterOptions,fileName_HCTSA,varargin)
% TS_Normalize  Trims and normalizes data from an hctsa analysis.
%
% Reads in data from HCTSA.mat, writes a trimmed, normalized version to HCTSA_N.mat
% For many normalization settings, each feature is normalized to the unit interval
% for the purposes of visualization and clustering.
%
% Time series with too few good values are removed (using each time series' own
% values). The features are then filtered and normalized (BF_FitTransform):
% by default, this is fitted on all time series; with 'trainIDs', it is fitted
% on those time series only and applied unchanged to all of them, so that the
% other time series can be used as an unseen test set.
%
%---USAGE:
% TS_Normalize('mixedSigmoid',[0.70,1]);
% TS_Normalize('mixedSigmoid',[0.70,1],'HCTSA.mat',false,false);
% TS_Normalize('mixedSigmoid',[0.70,1],'HCTSA.mat','trainIDs',trainIDs);
%
%---INPUTS:
% normFunction: String specifying how to normalize the data: 'none' (no
%               normalization), or a normalization method of BF_NormalizeMatrix
%               (e.g., 'mixedSigmoid', 'scaledRobustSigmoid', 'zscore', 'maxmin').
%               Default: 'mixedSigmoid'.
%
% filterOptions: Vector specifying thresholds for the minimum proportion of good
%                values required in a given row or column, in the form of a 2-vector:
%                [row proportion, column proportion]. If one of the filterOptions
%                is set to 1, will have no bad values in your matrix.
%                Default: [0.70,1].
%
% fileName_HCTSA: Custom hctsa data to import. Default: loaded from 'HCTSA.mat'.
%
% classVarFilter: whether to filter on zero variance of any given class (which
%                 can cause problems for many classification algorithms).
%                 Default: false.
%
% keepCalcTime: whether to keep TS_CalcTime. Default: false.
%
% (name-value option, after the inputs above):
% 'trainIDs': IDs of the time series to fit the feature filter and normalization
%             on (e.g., a training set). It is then applied to all time series.
%             Default: [] (all time series).
%
%---OUTPUT:
% outputFileName, the file the trimmed, normalized data is saved to (the input
%       filename with _N appended). Its normalizationInfo structure records the
%       settings (normFunction, filterOptions, classVarFilter, codeToRun) and, in
%       featureParams, a table with a row for each feature kept (ID and Name,
%       matching Operations) with the parameters of its normalizing
%       transformation (cf. the normParams output of BF_NormalizeMatrix).
%       These parameters are used by TS_Predict to normalize new data in the
%       same way as the training data. trainIDs lists the IDs of the time series
%       the filter and parameters were fitted on (all time series kept, by
%       default), and transform is the fitted transformation (the output of
%       BF_FitTransform, which can be applied to other data by BF_ApplyTransform).

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% --------------------------------------------------------------------------
%% Check Inputs
% --------------------------------------------------------------------------
if nargin < 1 || isempty(normFunction)
    fprintf(1,'Using the mixed sigmoidal transform: ''mixedSigmoid''\n')
    normFunction = 'mixedSigmoid';
end

if nargin < 2 || isempty(filterOptions)
    filterOptions = [0.70,1];
    % By default remove less than 70%-good-valued time series, & then less than
    % 100%-good-valued operations.
end
if any(filterOptions > 1)
    error('Set filterOptions as a length-2 vector with elements in the unit interval');
end
fprintf(1,['Removing time series with more than %.2f%% special-valued outputs\n' ...
            'Removing operations with more than %.2f%% special-valued outputs\n'], ...
            (1-filterOptions(1))*100,(1-filterOptions(2))*100);

% By default, work on HCTSA.mat, e.g., generated by TS_Init
if nargin < 3 || isempty(fileName_HCTSA)
    fileName_HCTSA = 'HCTSA.mat';
end

% Optional inputs classVarFilter and keepCalcTime, followed by name-value options:
isName = cellfun(@(x) ischar(x) || isstring(x),varargin);
firstName = find(isName,1);
if isempty(firstName)
    firstName = length(varargin) + 1;
end
if firstName > 3
    error('Too many inputs before the name-value options');
end
if firstName > 1 && ~isempty(varargin{1})
    classVarFilter = varargin{1};
else
    classVarFilter = false; % don't filter on individual class variance > 0 by default
end
if firstName > 2 && ~isempty(varargin{2})
    keepCalcTime = varargin{2};
else
    % Save space in normalized dataset
    keepCalcTime = false;
end
inputP = inputParser;
addParameter(inputP,'trainIDs',[],@(x) isempty(x) || isnumeric(x));
parse(inputP,varargin{firstName:end});
trainIDs = inputP.Results.trainIDs;

% --------------------------------------------------------------------------
%% Read data from local files
% --------------------------------------------------------------------------

% Load data:
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(fileName_HCTSA);
TS_Quality = TS_GetFromData(fileName_HCTSA,'TS_Quality');
MasterOperations = TS_GetFromData(fileName_HCTSA,'MasterOperations');

if keepCalcTime
    TS_CalcTime = TS_GetFromData(fileName_HCTSA,'TS_CalcTime');
end

% Maybe we kept the git repository info
gitInfo = TS_GetFromData(fileName_HCTSA,'gitInfo');

%-------------------------------------------------------------------------------
% In this script, each of these pieces of data (from the database) will be
% trimmed and normalized, and then saved to HCTSA_N.mat
%-------------------------------------------------------------------------------

% --------------------------------------------------------------------------
%% Trim down bad rows/columns
% --------------------------------------------------------------------------

% (i) NaNs in TS_DataMat mean values uncalculated in the matrix.
TS_DataMat(~isfinite(TS_DataMat)) = NaN; % Convert all nonfinite values to NaNs for consistency
% Need to also incorporate knowledge of bad entries in TS_Quality and filter these out:
TS_DataMat(TS_Quality > 0) = NaN;
numSpecialValues = nnz(TS_Quality > 0);
fprintf(1,'\nThere are %u special values in the data matrix.\n',numSpecialValues);
percGoodRows = mean(~isnan(TS_DataMat),2)*100;
percGoodCols = mean(~isnan(TS_DataMat),1)*100;
fprintf(1,'(pre-filtering): Time series vary from %.2f--%.2f%% good values\n',...
                min(percGoodRows),max(percGoodRows));
fprintf(1,'(pre-filtering): Features vary from %.2f--%.2f%% good values\n',...
                min(percGoodCols),max(percGoodCols));

% Now that all bad values are NaNs, and we can get on with the job of filtering them out

% (*) Filter based on proportion of bad entries. If either threshold is 1,
% the resulting matrix is guaranteed to be free from bad values entirely.

% Filter time series (rows)
keepRows = filterNaNs(TS_DataMat,filterOptions(1),'time series');
if any(~keepRows)
    fprintf(1,'Time series removed: %s.\n\n',BF_cat(TimeSeries.Name(~keepRows),','));
    TS_DataMat = TS_DataMat(keepRows,:);
    TS_Quality = TS_Quality(keepRows,:);
    TimeSeries = TimeSeries(keepRows,:);
    if keepCalcTime
        TS_CalcTime = TS_CalcTime(keepRows,:);
    end
end

%-------------------------------------------------------------------------------
%% Time series used to fit the feature filter and normalization
%-------------------------------------------------------------------------------
if height(TimeSeries)==1
    % When there is only a single time series, it doesn't actually make sense to normalize
    error('Only a single time series remains in the dataset -- normalization cannot be applied');
end
if isempty(trainIDs)
    isTrain = true(height(TimeSeries),1);
else
    isTrain = ismember(TimeSeries.ID,trainIDs);
    isRemoved = ~ismember(trainIDs,TimeSeries.ID);
    if any(isRemoved)
        warning('%u of the %u trainIDs are not among the time series kept (not in the data, or removed by filtering)',...
                        sum(isRemoved),length(trainIDs));
    end
    fprintf(1,['Fitting the feature filter and normalization on %u training time series' ...
                ' (and applying it to all %u time series)\n'],sum(isTrain),height(TimeSeries));
end

% --------------------------------------------------------------------------
%% Filter and normalize the features (operations)
% --------------------------------------------------------------------------
% Removes features with too few good values, constant features, (if
% classVarFilter) features constant within a class, and features that are all-NaN
% or constant after normalization; all fitted on the training time series only.
if classVarFilter && ~ismember('Group',TimeSeries.Properties.VariableNames)
    error('Group labels not assigned to time series, so cannot filter on class variance');
end
normSpec = struct('normFunction',normFunction,'filterOptions',filterOptions,...
                        'classVarFilter',classVarFilter);
if classVarFilter
    yTrain = TimeSeries.Group(isTrain);
else
    yTrain = [];
end
transform = BF_FitTransform(TS_DataMat(isTrain,:),normSpec,Operations,yTrain,true);
TS_DataMat = BF_ApplyTransform(TS_DataMat,transform,Operations);

% Keep the features of the transformation:
[~,keepCols] = ismember(transform.featureParams.ID,Operations.ID);
TS_Quality = TS_Quality(:,keepCols);
Operations = Operations(keepCols,:);
if keepCalcTime
    TS_CalcTime = TS_CalcTime(:,keepCols);
end

numBadEntries = sum(isnan(TS_DataMat(:)));
percBadEntries = numBadEntries/length(TS_DataMat(:))*100;
if numBadEntries==0
    fprintf(1,'No special-valued entries in the %ux%u data matrix!\n', ...
                size(TS_DataMat,1),size(TS_DataMat,2));
else
    fprintf(1,'%u special-valued entries (%4.2f%%) in the %ux%u data matrix.\n', ...
                numBadEntries,percBadEntries,size(TS_DataMat,1),size(TS_DataMat,2));
end

% ------------------------------------------------------------------------------
% Set default clustering details
% ------------------------------------------------------------------------------
ts_clust = struct('distanceMetric','none','Dij',[],...
                'ord',1:size(TS_DataMat,1),'linkageMethod','none');
op_clust = struct('distanceMetric','none','Dij',[],...
                'ord',1:size(TS_DataMat,2),'linkageMethod','none');

% --------------------------------------------------------------------------
%% Save results to file
% --------------------------------------------------------------------------
% Make a structure with statistics on normalization:
% Save the codeToRun, so you can check the settings used to run the normalization
% At the moment, only saves the first two arguments
codeToRun = sprintf('TS_Normalize(''%s'',[%f,%f])',normFunction, ...
                                        filterOptions(1),filterOptions(2));
% Store the features kept, with the parameters of their normalizing transformation,
% and the time series these were fitted on:
normalizationInfo = struct('normFunction',normFunction,'filterOptions', ...
                                    filterOptions,'codeToRun',codeToRun,...
                                    'featureParams',transform.featureParams,...
                                    'classVarFilter',classVarFilter,...
                                    'trainIDs',TimeSeries.ID(isTrain),...
                                    'transform',transform);

outputFileName = [fileName_HCTSA(1:end-4),'_N.mat'];

fprintf(1,'Saving the trimmed, normalized data to %s...',outputFileName);
if keepCalcTime
    save(outputFileName,'TS_DataMat','TS_Quality','TS_CalcTime','TimeSeries',...
            'Operations','MasterOperations','normalizationInfo',...
            'gitInfo','ts_clust','op_clust','-v7.3');
else
    save(outputFileName,'TS_DataMat','TS_Quality','TimeSeries',...
            'Operations','MasterOperations','normalizationInfo',...
            'gitInfo','ts_clust','op_clust','-v7.3');
end
fprintf(1,' Done.\n');

% Check whether output to screen is required:
if nargout == 0
    clear('outputFileName');
end

%-------------------------------------------------------------------------------
function keepInd = filterNaNs(XMat,nan_thresh,objectName)
    % Returns an index of rows of XMat with at least nan_thresh good values.

    if nan_thresh == 0
        keepInd = true(size(XMat,1),1);
        return
    else
        propNaN = mean(isnan(XMat),2); % proportion of NaNs across rows
        keepInd = (1-propNaN >= nan_thresh);
        if all(~keepInd)
            error('No %s had more than %4.2f%% good values.\nSet a more lenient threshold.',...
                                objectName,nan_thresh*100)
        end
        if all(keepInd)
            fprintf(1,['All %u %s have greater than %4.2f%% good values.' ...
                            ' Keeping them all.\n'], ...
                            length(keepInd),objectName,nan_thresh*100);
        else
            fprintf(1,['Removing %u %s with fewer than %4.2f%% good values:'...
                        ' from %u to %u.\n'],sum(~keepInd),objectName,...
                        nan_thresh*100,length(keepInd),sum(keepInd));
        end
    end
end
%-------------------------------------------------------------------------------

end
