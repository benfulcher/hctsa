function TS_plot_pca(whatData,TsorOps,showDist,classMeth,annotateParams)
% TS_plot_pca   2-dimensional feature-based representation of a time-series dataset.
%
% The low-dimensional representation is computed using PCA.
%
%---EXAMPLE USAGE:
%
% TSQ_plot_pca('norm');

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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
    fprintf(1,'Getting data from HCTSA_N\n');
end

if nargin < 2 || isempty(TsorOps)
    fprintf(1,'Using for time series\n');
    TsorOps = 'ts';
end
if ~any(ismember(TsorOps,{'ops','ts'}));
    error('Specify either operations (''ops'') or time series (''ts'').');
end

if nargin < 3 || isempty(showDist)
    showDist = 1;
end

if nargin < 4 || isempty(classMeth)
    classMeth = 'linclass';
end

if nargin < 5 || isempty(annotateParams)
    % Annotate 10 points by default
    fprintf(1,'Annotating 6 points by default with time series segment and names\n');
    annotateParams = struct('n',6,'textAnnotation','Name');
end
if ~isstruct(annotateParams)
    annotateParams = struct('n',annotateParams);
end

% ------------------------------------------------------------------------------
%% Load the data and group labeling from file
% ------------------------------------------------------------------------------
if strcmp(whatData,'cl') || strcmp(whatData,'norm')  || ischar(whatData)
    % Retrive data from local files
    if ismember(whatData,{'norm','cl'})
        whatDataFile = 'HCTSA_N.mat';
    else
        % Provided a custom filename to a datafile
        whatDataFile = whatData;
    end

    % Load in data:
    [TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

    % Construct labels for data points:
    if strcmp(TsorOps,'ts')
        dimensionLabels = {Operations.Name}; clear Operations % We just need their names
        if isstruct(annotateParams) || length(annotateParams) > 1 || annotateParams > 0
            dataLabels = {TimeSeries.Name};
            data_ids = [TimeSeries.ID];
            timeSeriesData = {TimeSeries.Data};
            if isfield(TimeSeries,'Group')
                load(whatDataFile,'groupNames')
                dataGroups = [TimeSeries.Group];
            else
                % No groups assigned:
                dataGroups = {};
                groupNames = {'all data'};
                % fprintf(1,'\n');
                % error('No groups assigned -- Use TS_LabelGroups to provide group information.')
            end
        end
    else
        dimensionLabels = {TimeSeries.Name}; clear TimeSeries
    end
    clear('TimeSeries','Operations'); % we no longer need you
else
    % The user provided data themself
    if ~isfield(whatData,'DataMat') && ~isfield(whatData,'TS_DataMat')
        error('No field ''DataMat'' (or ''TS_DataMat'') provided in the data input')
    elseif ~isfield(whatData,'Groups')
        error('No field ''Groups'' provided in the data input')
    end
    if isfield(whatData,'DataMat')
        TS_DataMat = whatData.DataMat;
    else
        TS_DataMat = whatData.TS_DataMat;
    end
    dataGroups = whatData.Groups;
    if isfield(whatData,'DimLabels')
        dimensionLabels = whatData.DimLabels;
    else
        dimensionLabels = {};
    end
    if isfield(whatData,'groupNames')
        groupNames = whatData.groupNames;
    else
        groupNames = {};
    end
    if isfield(whatData,'dataLabels')
        dataLabels = whatData.dataLabels;
    else
        dataLabels = {};
    end
    if isfield(whatData,'timeSeriesData')
        timeSeriesData = whatData.timeSeriesData;
    else
        timeSeriesData = {};
    end
end

if strcmp(TsorOps,'ops')
    % Take the transpose of the input data matrix for operations
    TS_DataMat = TS_DataMat';
end


% Label groups
if ~isempty(dataGroups)
    groupIndices = BF_ToGroup(dataGroups);
else
    groupIndices = {1:size(TS_DataMat,1)};
end
numGroups = length(groupIndices); % Number of groups

% ------------------------------------------------------------------------------
%% Do the dimensionality reduction using Matlab's built-in PCA algorithm
% ------------------------------------------------------------------------------
% Can't run PCA on data containing NaNs:
if any(isnan(TS_DataMat(:)))
    error('Data matrix contains NaNs');
end
fprintf(1,'Calculating principal components of the %u x %u data matrix...', ...
                    size(TS_DataMat,1),size(TS_DataMat,2));

% The new pca function is a much faster implementation to compute just the first
% 2 components:
[pcCoeff, pcScore, latent, ~, percVar] = pca(zscore(TS_DataMat),'NumComponents',2);
fprintf(1,' Done.\n');

% Work out the main contributions to each principle component
featLabel = cell(2,2); % Feature label :: proportion of total (columns are PCs)

% ------------------------------------------------------------------------------
% Plot this two-dimensional representation of the data using TS_plot_2d
% ------------------------------------------------------------------------------

nameString = 'PC';
for i = 1:2
    dataInfo.labels{i} = sprintf('%s %u (%.2f%%)',nameString,i,percVar(i));
end

dataInfo.GroupNames = groupNames;
dataInfo.GroupIndices = groupIndices;
dataInfo.DataLabels = dataLabels;
dataInfo.TimeSeriesData = timeSeriesData;

TS_plot_2d(pcScore(:,1:2),dataInfo,{},annotateParams,showDist,classMeth);

end
