function TS_plot_pca(whatData,showDist,classMeth,annotateParams)
% TS_plot_pca   2-dimensional feature-based representation of a time-series dataset.
%
% The low-dimensional representation is computed using PCA.
%
%---INPUTS:
% whatData, the hctsa data file (or structure) to use (input to TS_LoadData)
% showDist, whether to also plot marginal class distributions of each PC
% classMeth, the classification method when running classification in the PC space
% annotateParams, the annotation parameters used when plotting the dataset using
%                 TS_plot_2d. Can also specify an integer; will annotate this
%                 many time series segments to the 2d scatter plot.
%
%---EXAMPLE USAGE:
% (*) Plot a PCA projection of the normalized data stored in HCTSA_N.mat
% >> TSQ_plot_pca('norm');

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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm'; % normalized data by default, from HCTSA_N.mat
end

if nargin < 2 || isempty(showDist)
    showDist = 1;
end

if nargin < 3 || isempty(classMeth)
    classMeth = 'linclass';
end

if nargin < 4 || isempty(annotateParams)
    % Annotate 6 points by default
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
[TS_DataMat,TimeSeries,~,theFile] = TS_LoadData(whatData);

% Retrieve group names also:
fileVarsStruct = whos('-file',theFile);
fileVars = {fileVarsStruct.name};
if ismember('groupNames',fileVars)
    groupNames = load(theFile,'groupNames');
    groupNames = groupNames.groupNames;
else
    groupNames = {};
end

% ------------------------------------------------------------------------------
%% Do the dimensionality reduction using Matlab's built-in PCA algorithm
% ------------------------------------------------------------------------------
fprintf(1,'Calculating principal components of the %u x %u data matrix...\n', ...
                    size(TS_DataMat,1),size(TS_DataMat,2));

% Use pca to compute the first two principal components:
if ~any(isnan(TS_DataMat))
    [~,pcScore,~,~,percVar] = pca(TS_DataMat,'NumComponents',2,'Centered','on');
else
    warning(sprintf(['Data matrix contains %.2g%% NaNs. Attempting to estimate covariances on remaining data...\n' ...
                '(Could take some time...)'],100*mean(isnan(TS_DataMat(:)))))
    % Data matrix contains NaNs; try the pairwise rows approximation to the
    % covariance matrix:
    [~,pcScore,~,~,percVar] = pca(TS_DataMat,'Rows','pairwise','Centered','on');
    % If this fails (covariance matrix not positive definite), could try
    % the 'algorithm','als' option in pca...
end
fprintf(1,'---Done.\n');

% ------------------------------------------------------------------------------
% Plot this two-dimensional representation of the data using TS_plot_2d
% ------------------------------------------------------------------------------

% Set feature labels:
nameString = 'PC';
featureLabels = cell(2,1);
for i = 1:2
    featureLabels{i} = sprintf('%s %u (%.2f%% var)',nameString,i,percVar(i));
end

TS_plot_2d(pcScore(:,1:2),TimeSeries,featureLabels,groupNames,annotateParams,showDist,classMeth)

end
