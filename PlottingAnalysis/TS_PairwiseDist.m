function Dij = TS_PairwiseDist(tsOrOps,whatData,distanceMetric,doSave)
% TS_PairwiseDist   Pairwise distances between either TimeSeries or Operations
%
% (If doSave), stores pairwise distances back the data back in the local .mat
% file. Warning: this can be huge for large datasets.
%
%---INPUTS:
%
% tsOrOps: Compute pairwise distances between all pairs of time series
%               ('ts', default), or operations ('ops')
%
% whatData: what hctsa data source to use (cf. TS_LoadData)
%
% distanceMetric: what distance metric to use, e.g., 'euclidean' (default for
%                 time series) or 'corr_fast' (default for operations)

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs, set defaults:
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(tsOrOps)
    tsOrOps = 'ts';
    fprintf(1,'Time series by default.\n');
end
if ~ismember(tsOrOps,{'ts','ops'})
    error('Invalid specifier ''%s'' (should be ''ts'' or ''ops'')',tsOrOps);
end

if nargin < 2
    whatData = 'norm';
end

if nargin < 3
    switch tsOrOps
    case 'ts'
        distanceMetric = 'euclidean';
    case 'ops'
        distanceMetric = 'corr_fast';
    end
end

% doSave
if nargin < 4
    % By default, save back to the file specified as whatData
    doSave = 1;
end

% ------------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------------

if ischar(whatData)
    [TS_DataMat,~,~,theFile] = TS_LoadData(whatData);
else
    % If provided a matrix, assume it's from HCTSA_N.mat by default:
    TS_DataMat = whatData;
    theFile = 'HCTSA_N.mat';
end

% ------------------------------------------------------------------------------
% Compute pairwise distances between all objects
% ------------------------------------------------------------------------------
if strcmp(tsOrOps,'ops')
    % Transpose if we want columns ('ops')
    TS_DataMat = TS_DataMat';
end

fprintf(1,'Computing %s distances between %u objects...',distanceMetric,size(TS_DataMat,1));
Dij = BF_pdist(TS_DataMat,distanceMetric,1,[],1);
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% Save back to the local hctsa .mat file
% ------------------------------------------------------------------------------
if doSave
    switch tsOrOps
    case 'ts'
        theWhat = 'time series';
    case 'ops'
        theWhat = 'operations';
    end
    fprintf(1,'Saving pairwise distance information for %s back to %s...',...
                        theWhat,theFile);

    switch tsOrOps
    case 'ts'
        ts_clust = struct('Dij',Dij,'distanceMetric',distanceMetric);
        save(theFile,'ts_clust','-append')
    case 'ops'
        op_clust = struct('Dij',Dij,'distanceMetric',distanceMetric);
        save(theFile,'op_clust','-append')
    end
    
    fprintf(1,' Done.\n');
end

end
