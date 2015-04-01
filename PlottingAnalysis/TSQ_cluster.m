% --------------------------------------------------------------------------
% TSQ_cluster
% --------------------------------------------------------------------------
% 
% Reads in normalized data from HCTSA_N.mat, clusters the data matrix by
% reordering rows and columns with linkage clustering, and then saves the result
% as HCTSA_cl.mat
% 
%---EXAMPLE USAGE:
% TSQ_cluster;
% 
%---INPUTS:
% distanceMetricRow: specifies the distance metric for computing distances
%                   between time series (rows)
% linkageMethodRow: specifies the linkage method for clustering time series (rows)
% distanceMetricCol: distance metric for distances between operations (columns)
% linkageMethodCol: linkage method for clustering operations (columns)
% subSet: should be a cell with three components: the first should be either
%         'cl' or 'norm' for the row and column indicies provided in the 2nd and
%         3rd components of the cell. e.g., {'norm',[1,5,...],[]} will load in
%         HCTSA_N, and cluster the subset of rows [1,5,...] and all columns.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TSQ_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol, subSet)

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

% use same method as for rows
if nargin < 3 || isempty(distanceMetricCol)
    % fprintf(1,'Clustering columns in the same way as for the rows: using %s\n',clusterMethRow)
    distanceMetricCol = 'corr_fast';
end

% clustering settings for columns
if nargin < 4 || isempty(linkageMethodCol)
    linkageMethodCol = 'average';
end

% Subsets -- only cluster a subset of the full data matrix
if nargin < 5
    subSet = []; % cluster the full input matrix by default
end

if ~isempty(subSet)
    if (length(subSet) ~= 3)
        error('The subset should specify ''norm'' or ''cl'' and the subset.')
    elseif ~ismember(subSet{1},{'norm','cl'})
        error('The first component of subset should be either ''norm'' or ''cl''.')
    end
end

% --------------------------------------------------------------------------
% Record information about clustering settings so can run the same back again if needed later:
% --------------------------------------------------------------------------
codeToCluster = ['TSQ_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, ' ...
                                    'linkageMethodCol, subSet)'];
clusteringInfo = struct('distanceMetricRow',distanceMetricRow,'linkageMethodRow', ...
                linkageMethodRow,'distanceMetricCol',distanceMetricCol,'linkageMethodCol',linkageMethodCol, ...
                'Subset',subSet,'codeToCluster',codeToCluster);

% --------------------------------------------------------------------------
%% Load information from local files
% --------------------------------------------------------------------------
if isempty(subSet) || strcmp(subSet{1},'norm')
    TheFile = 'HCTSA_N.mat';
else % subset of the clustered matrix
    TheFile = 'HCTSA_cl.mat';
end
wn = which(TheFile); % Check that HCTSA_N exists
if isempty(wn);
    error('%s not found.',TheFile);
end
fprintf(1,'Loading data from %s...',TheFile)
load('HCTSA_N.mat','TS_DataMat','TimeSeries','Operations','MasterOperations','normalizationInfo')
fprintf(1,' Done.\n')

% --------------------------------------------------------------------------
%% Take subsets of the data
% --------------------------------------------------------------------------
if ~isempty(subSet)
    fprintf(1,'We are now implementing subset behaviour...\n')
    % subs is in the form {[rowrange],[columnrange]}; a cell of two vectors

    if ~isempty(subSet{2}); % row subset
       r = subSet{2};
       TS_DataMat = TS_DataMat(r,:);
       TimeSeries = TimeSeries(r);
       fprintf(1,'Reduced rows/timeseries from %u to %u according to the specified subset.\n',length(TimeSeries),length(r));
    end

    if ~isempty(subSet{3}); % column subset
        r = subSet{3};
        TS_DataMat = TS_DataMat(:,r);
        Operations = Operations(r);
        fprintf(1,'Reduced columns/operations from %u to %u according to the specified subset\n',length(Operations),length(r));
    end
end

% --------------------------------------------------------------------------
%% Do the clustering
% --------------------------------------------------------------------------
fprintf(1,'Clustering the full %u x %u data matrix.\n',length(TimeSeries),length(Operations))

% Cluster rows
if ~(ischar(distanceMetricRow) && ismember(distanceMetricRow,{'none','nothing'})) && size(TS_DataMat,1) > 1
    % (can specify 'none' to do no clustering)
    fprintf(1,'\n----Clustering rows using a distance metric %s and %s linkage method----\n',...
                    distanceMetricRow,linkageMethodRow); tic
    % [~, ord_row] = TSQ_us_cluster(TS_DataMat,clusterMethRow,clusterParamsRow);
    ord_row = TSQ_ClusterReorder(TS_DataMat,distanceMetricRow,linkageMethodRow);
    fprintf(1,'Row clustering took %s.\n',BF_thetime(toc))
else
    ord_row = 1:size(TS_DataMat,1); % don't reorder at all
end

% Cluster columns
if ~(ischar(distanceMetricCol) && ismember(distanceMetricCol,{'none','nothing'})) && size(TS_DataMat,2) > 1
    % (can specify 'none' to do no clustering)
    fprintf(1,'\n----Clustering columns using a distance metric %s and %s linkage method----\n',...
                    distanceMetricCol,linkageMethodCol); tic
    ord_col = TSQ_ClusterReorder(TS_DataMat',distanceMetricCol,linkageMethodCol);
    % [~, ord_col] = TSQ_us_cluster(TS_DataMat',clusterMethCol,clusterParamsCol);
    fprintf(1,'Column clustering took %s.\n',BF_thetime(toc))
else
    ord_col = 1:size(TS_DataMat,2); % don't reorder at all
end

% Save to the clusteringInfo structure:
clusteringInfo.ord_row = ord_row;
clusteringInfo.ord_col = ord_col;

% --------------------------------------------------------------------------
%% Reorder data structures according to the clustering
% --------------------------------------------------------------------------

% Reorder data matrix
TS_DataMat = TS_DataMat(ord_row,ord_col);

% Reorder time series metadata
TimeSeries = TimeSeries(ord_row);

% Reorder operation metadata
Operations = Operations(ord_col);

% --------------------------------------------------------------------------
%% Save results to file
% --------------------------------------------------------------------------

% Save information about the clustering to the file
% You can run codeToCluster after loading all the clustering info using an eval statement...?
fprintf(1,'Saving the clustered data as ''HCTSA_cl''...')
save('HCTSA_cl.mat','TS_DataMat','TimeSeries','Operations','MasterOperations','normalizationInfo','clusteringInfo')
fprintf(1,' Done.\n');


end