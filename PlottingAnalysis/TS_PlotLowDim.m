function f = TS_PlotLowDim(whatData,whatAlgorithm,showDist,annotateParams,cfnParams)
% TS_PlotLowDim   2-dim feature-based representation of a time-series dataset
%
% The low-dimensional representation is computed using PCA by default.
%
%---INPUTS:
% whatData, the hctsa data file (or structure) to use (input to TS_LoadData)
% whatAlgorithm, what dimensionality-reduction algorithm to use ('pca' is default),
%                   also 'tSNE'.
% showDist, whether to also plot marginal class distributions of each PC
% annotateParams, the annotation parameters used when plotting the dataset using
%                 TS_Plot2d. Can also specify an integer; will annotate this
%                 many time series segments to the 2d scatter plot.
% cfnParams, parameters for classification in 2-d space (cf., GiveMeDefaultClassificationParams)
%
%---OUTPUT:
% Figure handle, f
%
%---EXAMPLE USAGE:
% (*) Plot a PCA projection of the normalized data stored in HCTSA_N.mat
% >> TS_PlotLowDim('norm','pca');

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm'; % load in normalized data by default, from HCTSA_N.mat
end
if nargin < 2 || isempty(whatAlgorithm)
    whatAlgorithm = 'pca';
end
if nargin < 3 || isempty(showDist)
    showDist = true;
end
if nargin < 4 || isempty(annotateParams)
    % Annotate 6 points by default
    fprintf(1,'Annotating 6 points by default with time series segment and names\n');
    annotateParams = struct('n',6,'textAnnotation','Name','userInput',false);
end
if ~isstruct(annotateParams)
    annotateParams = struct('n',annotateParams);
end

%-------------------------------------------------------------------------------
%% Load the data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);
numFeatures = height(Operations);

% (for classification in the embedding space):
if nargin < 5 || isempty(cfnParams)
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries);
end

% Assign group labels (removing unlabeled data):
[TS_DataMat,TimeSeries] = FilterLabeledTimeSeries(TS_DataMat,TimeSeries);
% Give basic info about the represented classes:
TellMeAboutLabeling(TimeSeries);

% Filter down a reduced feature set if required:
[TS_DataMat,Operations] = FilterFeatures(TS_DataMat,Operations,cfnParams);

% ------------------------------------------------------------------------------
% Do the dimensionality reduction
% ------------------------------------------------------------------------------
[lowDimComponents,componentLabels] = BF_LowDimProject(TS_DataMat,whatAlgorithm,2);

%-------------------------------------------------------------------------------
% Plot two-dimensional representation of the data using TS_Plot2d:
f = TS_Plot2d(lowDimComponents,TimeSeries,componentLabels,annotateParams,showDist,cfnParams);

%-------------------------------------------------------------------------------
% Clear output
if nargout==0
    clear('f')
end

end
