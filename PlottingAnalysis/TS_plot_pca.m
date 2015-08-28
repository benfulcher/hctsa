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

% Load in data:
[TS_DataMat,TimeSeries,Operations,theFile] = TS_LoadData(whatData);

% Retrieve group names also:
fileVarsStruct = whos('-file',theFile);
fileVars = {fileVarsStruct.name};
if ismember('groupNames',fileVars)
    groupNames = load(theFile,'groupNames');
    groupNames = groupNames.groupNames;
else
    groupNames = {};
end

if strcmp(TsorOps,'ops')
    % Take the transpose of the input data matrix for operations
    TS_DataMat = TS_DataMat';
end

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
featureLabels = cell(2,1);
for i = 1:2
    featureLabels{i} = sprintf('%s %u (%.2f%%)',nameString,i,percVar(i));
end

TS_plot_2d(pcScore(:,1:2),TimeSeries,featureLabels,groupNames,annotateParams,showDist,classMeth)

end
