function reducedIDs = TS_ReduceFeatureSet(whatData,distThreshold,useSpearman)
% Unsupervised clustering to a reduced feature set
%
%---Output:
% reducedIDs, A set of Operation IDs that forms the new reduced feature set

% ------------------------------------------------------------------------------
% Copyright (C) 2018, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1
    whatData = 'norm';
end
if nargin < 2
    distThreshold = 0.2;
end
if nargin < 3
    useSpearman = true;
end

%-------------------------------------------------------------------------------
% Load in data:
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

% Compute correlation distances between all pairs of features on the data:
if useSpearman
    rawDistVec = pdist(TS_DataMat','spearman');
else
    rawDistVec = pdist(TS_DataMat','corr');
end

% Cluster down with abscorr linkage clustering:
[distMat_cl,cluster_Groupi,ord,handles] = BF_ClusterDown(rawDistVec,...
                'clusterThreshold',distThreshold,...
                'whatDistance','abscorr',...
                'linkageMeth','complete');

% Determine cluster centers:
clusterCenterInd = cellfun(@(x)x(1),cluster_Groupi);

% Convert to feature IDs:
reducedIDs = Operations.ID(clusterCenterInd);

% Give user feedback:
fprintf(1,'Reduced to set of %u features using complete linkage clustering on abscorr distances (%.2f)\n',...
                length(reducedIDs),distThreshold);

end
