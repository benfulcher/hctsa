% --------------------------------------------------------------------------
% TSQ_cluster
% --------------------------------------------------------------------------
% 
% Reads in normalized data from HCTSA_N.mat, clusters the data matrix by
% reordering rows and columns, then saves the result as HCTSA_cl.mat
% 
%---INPUTS:
% ClusterMethRow: specifies the clustering method for rows/time series (default is 'linkage')
% ClusterParamsRow: specifies a cell of parameters specifying the clustering parameters,
%           including the distance metric, etc. (default is euclidean distances
%           and average linkage)
% ClusterMethCol: specifies the clustering method for columns/operations (default is the
%         same as rows)
% ClusterParamsCol: clustering settings for columns (default is correlation distances
%           and average linkage)
% SubSet: should be a cell with three components: the first should be either
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

function TSQ_cluster(ClusterMethRow, ClusterParamsRow, ClusterMethCol, ClusterParamsCol, SubSet)

% --------------------------------------------------------------------------
%% Check input arguments
% --------------------------------------------------------------------------
if nargin < 1 || isempty(ClusterMethRow)
    ClusterMethRow = 'linkage';
end

% Clustering settings for rows
% (can specify the string 'none' for no clustering)
if nargin < 2
    ClusterParamsRow = {'euclidean','average',0,[],0};
end

% use same method as for rows
if nargin < 3 || isempty(ClusterMethCol)
    fprintf(1,'Clustering columns in the same way as for the rows: using %s\n',ClusterMethRow)
    ClusterMethCol = ClusterMethRow;
end

% clustering settings for columns
if nargin < 4 || isempty(ClusterParamsCol)
    ClusterParamsCol = {'corr','average',0,[],0};
end

% Subsets -- only cluster a subset of the full data matrix
if nargin < 5
    SubSet = []; % cluster the full input matrix by default
end

if ~isempty(SubSet)
    if (length(SubSet) ~= 3)
        error('The subset should specify ''norm'' or ''cl'' and the subset.')
    elseif ~ismember(SubSet{1},{'norm','cl'})
        error('The first component of subset should be either ''norm'' or ''cl''.')
    end
end

% --------------------------------------------------------------------------
% Record information about clustering settings so can run the same back again if needed later:
% --------------------------------------------------------------------------
CodeToCluster = ['TSQ_cluster(ClusterMethRow, ClusterParamsRow, ClusterMethCol, ' ...
                                    'ClusterParamsCol, SubSet)'];
ClusteringInfo = struct('ClusterMethRow',ClusterMethRow,'ClusterParamsRow', ...
                ClusterParamsRow,'ClusterMethCol',ClusterMethCol,'ClusterParamsCol',ClusterParamsCol, ...
                'Subset',SubSet,'CodeToCluster',CodeToCluster);

% --------------------------------------------------------------------------
%% Load information from local files
% --------------------------------------------------------------------------
if isempty(SubSet) || strcmp(SubSet{1},'norm')
    TheFile = 'HCTSA_N.mat';
else % subset of the clustered matrix
    TheFile = 'HCTSA_cl.mat';
end
wn = which(TheFile); % Check that HCTSA_N exists
if isempty(wn);
    error('%s not found.',TheFile);
end
fprintf(1,'Loading data from %s...',TheFile)
load('HCTSA_N.mat','TS_DataMat','TimeSeries','Operations','MasterOperations','NormalizationInfo')
fprintf(1,' Done.\n')


% --------------------------------------------------------------------------
%% Take subsets of the data
% --------------------------------------------------------------------------
if ~isempty(SubSet)
    fprintf(1,'We are now implementing subset behaviour...\n')
    % subs is in the form {[rowrange],[columnrange]}; a cell of two vectors

    if ~isempty(SubSet{2}); % row subset
       r = SubSet{2};
       TS_DataMat = TS_DataMat(r,:);
       TimeSeries = TimeSeries(r);
       fprintf(1,'Reduced rows/timeseries from %u to %u according to the specified subset.\n',length(TimeSeries),length(r));
    end

    if ~isempty(SubSet{3}); % column subset
        r = SubSet{3};
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
if ~(ischar(ClusterMethCol) && ismember(ClusterMethRow,{'none','nothing'})) % can specify 'none' to do no clustering
    fprintf(1,'Clustering rows...\n'); tic
    [~, acgir] = TSQ_us_cluster(TS_DataMat,ClusterMethRow,ClusterParamsRow,'ts');
    fprintf(1,'Row clustering took %s.\n',BF_thetime(toc))
else
    acgir = {};
end

% Cluster columns
if ~(ischar(ClusterMethCol) && ismember(ClusterMethCol,{'none','nothing'})) && size(TS_DataMat,2)>1 % can specify 'none' to do no clustering
    fprintf(1,'Clustering columns...\n'); tic
    [~, acgic] = TSQ_us_cluster(TS_DataMat',ClusterMethCol,ClusterParamsCol,'mets');
    fprintf(1,'Column clustering took %s.\n',BF_thetime(toc))
else
    acgic = {};
end

% --------------------------------------------------------------------------
%% Reorder data structures according to the clustering
% --------------------------------------------------------------------------
% Get the permutation vectors ordr (reordering for rows) and ordc
% (reordering for columns)

if isempty(acgir)
    ordr = 1:size(TS_DataMat,1); % don't reorder at all
elseif iscell(acgir)
    ordr = vertcat(acgir{:});
else
    ordr = acgir;
end
if isempty(acgic);
    ordc = 1:size(TS_DataMat,2); % don't reorder at all
elseif iscell(acgic)
    ordc = vertcat(acgic{:});
else
    ordc = acgic;
end

% Reorder data matrix
TS_DataMat = TS_DataMat(ordr,ordc);

% Reorder time series metadata
TimeSeries = TimeSeries(ordr);

% Reorder operation metadata
Operations = Operations(ordc);


% --------------------------------------------------------------------------
%% Save results to file
% --------------------------------------------------------------------------

% Save information about the clustering to the file
% You can run CodeToCluster after loading all the clustering info using an eval statement...?
fprintf(1,'Saving the clustered data as ''HCTSA_cl''...')
save('HCTSA_cl.mat','TS_DataMat','TimeSeries','Operations','MasterOperations','NormalizationInfo','ClusteringInfo')
fprintf(1,' Done.\n');


end