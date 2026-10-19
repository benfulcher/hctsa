function reducedIDs = TS_ReduceFeatureSet(whatData,distThreshold,useSpearman)
% TS_ReduceFeatureSet   Unsupervised clustering to a reduced feature set
%
% Clusters features by the magnitude of their correlation across time series
% (complete-linkage clustering on the distance 1 - |rho|), and returns the
% most central feature of each cluster. Features that are constant (or all NaN)
% have no defined correlation and are left out. A figure showing the clustering
% is drawn (cf. BF_ClusterDown).
%
%---INPUTS (all optional; empty inputs take the defaults):
% whatData, the hctsa data to use (cf. TS_LoadData; default: 'norm')
% distThreshold, the distance, 1 - |rho|, at which the dendrogram is cut into clusters
%               (default: 0.2). All features in a cluster have |rho| >= 1 - distThreshold.
% useSpearman, true to use Spearman correlations (default), false for Pearson
%
%---OUTPUT:
% reducedIDs, A set of Operation IDs that forms the new reduced feature set
%             (the most central feature of each cluster, largest clusters first)

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

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
end
if nargin < 2 || isempty(distThreshold)
    distThreshold = 0.2;
end
if nargin < 3 || isempty(useSpearman)
    useSpearman = true;
end

%-------------------------------------------------------------------------------
% Load in data:
[TS_DataMat,~,Operations] = TS_LoadData(whatData);

% Leave out features with no defined correlation (constant or all NaN):
featureStd = std(TS_DataMat,0,1,'omitnan');
isDefined = (featureStd > 0);
if ~all(isDefined)
    fprintf(1,'Leaving out %u constant or all-NaN features.\n',sum(~isDefined));
end
Operations = Operations(isDefined,:);
TS_DataMat = TS_DataMat(:,isDefined);

% Compute correlation distances (1 - rho) between all pairs of features, using the
% overlapping good values of each pair when there are NaNs:
if useSpearman
    rawDistVec = BF_pdist(TS_DataMat','Spearman',true,[],true);
else
    rawDistVec = BF_pdist(TS_DataMat','Pearson',true,[],true);
end

% Cluster down with abscorr linkage clustering:
[~,cluster_Groupi] = BF_ClusterDown(rawDistVec,...
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
