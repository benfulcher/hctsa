function TS_cluster(distanceMetricRow,linkageMethodRow,distanceMetricCol,linkageMethodCol,doSave,theFile)
% TS_cluster    Linkage clustering for hctsa data.
%
% Reads in normalized data from HCTSA_N.mat, clusters the data matrix by
% reordering rows and columns with linkage clustering, and then saves the result
% back to HCTSA_N.mat (or a custom hctsa file).
%
%---EXAMPLE USAGE:
% TS_cluster;
%
%---INPUTS:
% distanceMetricRow: specifies the distance metric for computing distances
%                   between time series (rows)
% linkageMethodRow: specifies the linkage method for clustering time series (rows)
% distanceMetricCol: distance metric for distances between operations (columns)
% linkageMethodCol: linkage method for clustering operations (columns)
% doSave: controls whether to save pairwise distance information (which can be large)
%         to the local file (HCTSA_N.mat)

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

% --------------------------------------------------------------------------
%% Check input arguments:
% --------------------------------------------------------------------------
if nargin < 1 || isempty(distanceMetricRow)
    distanceMetricRow = 'euclidean';
end

% Clustering settings for rows
% (can specify the string 'none' for no clustering)
if nargin < 2 || isempty(linkageMethodRow)
    linkageMethodRow = 'average';
end

% Use fast approximations to correlations;
if nargin < 3 || isempty(distanceMetricCol)
    % fprintf(1,'Clustering columns in the same way as for the rows: using %s\n',clusterMethRow)
    distanceMetricCol = 'corr_fast';
end

% clustering settings for columns
if nargin < 4 || isempty(linkageMethodCol)
    linkageMethodCol = 'average';
end

% argument 5, doSave, set defaults later
if nargin >= 5
    if length(doSave)~=2
        error('doSave has two components: for (i) time series and (ii) operations');
    end
end

% What file to use:
if nargin < 6
    theFile = 'norm';
end

% --------------------------------------------------------------------------
%% Load information from local mat file
% --------------------------------------------------------------------------
% Interpret 'norm':
if strcmp(theFile,'norm')
    theFile = 'HCTSA_N.mat';
end

% Check that the file exists:
if ~exist(theFile,'file');
    error('%s not found.',theFile);
end
fileVarsStruct = whos('-file',theFile);
fileVars = {fileVarsStruct.name};
matrixSize = fileVarsStruct(strcmp(fileVars,'TS_DataMat')).size;
numTs = matrixSize(1);
numOps = matrixSize(2);

% ------------------------------------------------------------------------------
% doSave, whether to save full pairwise distance information to the local file
% for [ts,ops]
% ------------------------------------------------------------------------------
% By default, store full pairwise distance information only if fewer than 1000 objects
if nargin < 5 || isempty(doSave)
    doSave = [1,1];
    if numTs > 1000
        doSave(1) = 0;
    end
    if numOps > 1000
        doSave(2) = 0;
    end
end

% ------------------------------------------------------------------------------
% Compute and store pairwise distances, then cluster:
% ------------------------------------------------------------------------------
fprintf(1,'Clustering the full %u x %u data matrix.\n',numTs,numOps);

% Time series:
[ts_ord,ts_update] = getDistances(distanceMetricRow,'ts',linkageMethodRow,numTs,doSave(1));

% Operations:
[op_ord,op_update] = getDistances(distanceMetricCol,'ops',linkageMethodCol,numOps,doSave(2));

% --------------------------------------------------------------------------
%% Save results to file
% --------------------------------------------------------------------------
% Save information about the clustering to the file
% You can run codeToCluster after loading all the clustering info using an eval statement...?
if ~ts_update && ~op_update
    fprintf(1,['Data already clustered:\nTime series by %s / %s linkage\n' ...
                'Operations by %s / %s linkage.\n'],...
                distanceMetricRow,linkageMethodRow,distanceMetricCol,linkageMethodCol);
end

% Time series clustering needs to be updated
if ts_update
    fprintf(1,'Saving the clustering information for time series back to %s...',theFile);
    load(theFile,'ts_clust');
    ts_clust.ord = ts_ord;
    ts_clust.linkageMethod = linkageMethodRow;
    save(theFile,'-append','ts_clust');
    fprintf(1,' Done.\n');
end

% Operation clustering needs to be updated
if op_update
    fprintf(1,'Saving the clustering information for operations back to %s...',theFile);
    load(theFile,'op_clust');
    op_clust.ord = op_ord;
    op_clust.linkageMethod = linkageMethodCol;
    save(theFile,'-append','op_clust');
    fprintf(1,' Done.\n');
end


% ------------------------------------------------------------------------------
function [ordering,didUpdate] = getDistances(dMetricNow,tsOrOps,theLinkageMethod,numObjects,saveBack)

    switch tsOrOps
    case 'ts'
        theWhat = 'time series';
        theVar = 'ts_clust';
    case 'ops'
        theWhat = 'operations';
        theVar = 'op_clust';
    end

    % Set to zero if no update was performed for (1) distances, and (2) linkage clustering
    didUpdate = ones(2,1);

    if ~(ischar(dMetricNow) && ismember(dMetricNow,{'none','nothing'}))
        % (can specify 'none' to do no clustering)

        % Check if this already exists:
        if ismember(theVar,fileVars)
            % Clustering has already been performed:
            clust_data = load(theFile,theVar);
            clust_data = clust_data.(theVar);
            if strcmp(clust_data.distanceMetric,dMetricNow)
                fprintf(1,'Loaded %s distances (%s) stored in %s\n',...
                            theWhat,clust_data.distanceMetric,theFile);
                Dij = clust_data.Dij;
                didUpdate(1) = 0;
            else
                fprintf(1,'Changing distance metric for %s from ''%s'' to ''%s''.\n',...
                            theWhat,clust_data.distanceMetric,dMetricNow);
                Dij = TS_PairwiseDist(tsOrOps,theFile,dMetricNow,saveBack);
            end
        else
            Dij = TS_PairwiseDist(tsOrOps,theFile,dMetricNow,saveBack);
        end

        % Reorder using linkage clustering:
        if didUpdate(1)==0 && isfield(clust_data,'linkageMethod') && ...
                strcmp(clust_data.linkageMethod,theLinkageMethod)
            % Already been clusterd previously under the same distance metric and linkage method:
            ordering = clust_data.ord;
            didUpdate(2) = 0;
        else
            ordering = BF_ClusterReorder([],Dij,theLinkageMethod);
        end
    else
        ordering = 1:numObjects; % don't reorder at all
    end

    didUpdate = any(didUpdate);
end
% ------------------------------------------------------------------------------

end
