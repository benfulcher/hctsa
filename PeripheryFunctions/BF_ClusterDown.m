function [distMat_cl,cluster_Groupi,ord,handles] = BF_ClusterDown(distMat,varargin)
% BF_ClusterDown    Reduce a pairwise distance matrix into smaller clusters.
%
% Yields a visualization of a pairwise distance matrix, including a set of
% clusters of objects showing highly correlated behavior.
%
%---INPUTS:
% distMat, a pairwise distance vector or matrix (e.g., generated from pdist)
% [EXTRA OPTIONS]:
% 'clusterThreshold', the threshold in the dendrogram for forming clusters
% 'objectLabels', labels for axes
% 'whatDistance', the type of distance metric used:
%        (i) 'corr' (default): provide 1-R, where R is correlation coefficient
%        (ii) 'abscorr': provide correlation distances, 1-R, but clustering is
%                        done on abscorr distances (sign of R is ignored), and
%                        correlations, R, are plotted
%        (iii) 'abscorr_ii': provide correlation distances, 1-R, but clustering
%                        is done on abscorr distances(sign of R is ignored), and
%                        absolute correlations, |R|, are plotted
%        (iv) 'general': provide any distance matrix, plots the distance matrix
% 'linkageMeth', the linkage method to use (for Matlab's linkage function)
%
%---OUTPUTS:
% distMat_cl, the distance matrix re-ordered by clustering
% cluster_Groupi, the assignment of nodes to clusters (cell of indices)
% ord, the ordering of objects in the dendrogram
% handles, structure containing handles to plotting elements, f (figure),
%               ax1 (dendrogram), ax2 (pairwise correlation matrix)

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs
%-------------------------------------------------------------------------------
% Custom plotting options using an input parser
inputP = inputParser;

% Cluster threshold:
default_clusterThreshold = 0.2;
addParameter(inputP,'clusterThreshold',default_clusterThreshold,@isnumeric);

% Object labels for axes
default_objectLabels = {};
addParameter(inputP,'objectLabels',default_objectLabels,@iscell);

% Default input is a distance matrix based on correlation coefficients
default_whatDistance = 'corr';
addParameter(inputP,'whatDistance',default_whatDistance,@ischar);

% errth Set error threshold instead of just top number
% default_errTh = 30;
% addParameter(inputP,'errTh',default_errTh,@isnumeric);

% miOrCorr
% default_miOrCorr = 'corr'; % 'mi','corr'
% addParameter(inputP,'miOrCorr',default_miOrCorr,@ischar);

% clusterMeth, what method to use to do the clustering
% default_clusterMeth = 'linkage'; % 'linkage', 'kmedoids'
% addParameter(inputP,'clusterMeth',default_clusterMeth,@ischar);

% linkageMeth, what method to use to do the linkage clustering
default_linkageMeth = 'average';
addParameter(inputP,'linkageMeth',default_linkageMeth,@ischar);

% plotBig, whether to make a big plot
% default_plotBig = 0;
% check_plotBig = @(x) x==1 || x==0;
% addParameter(inputP,'plotBig',default_plotBig,check_plotBig);

%% Parse inputs:
parse(inputP,varargin{:});

% Make variables from results of input parser:
clusterThreshold = inputP.Results.clusterThreshold;
objectLabels = inputP.Results.objectLabels;
whatDistance = inputP.Results.whatDistance;
% errTh = inputP.Results.errTh;
% miOrCorr = inputP.Results.miOrCorr;
% clusterMeth = inputP.Results.clusterMeth;
linkageMeth = inputP.Results.linkageMeth;
% plotBig = inputP.Results.plotBig;
clear inputP;

%-------------------------------------------------------------------------------
% Work with vector form of distances (hopefully enough memory)
%-------------------------------------------------------------------------------
if any(size(distMat)==1)
    distVec = distMat;
else
    distVec = squareform(distMat); % convert for vector version
end
% Solve for number of items from quadratic formula:
numItems = (1+sqrt(1+8*length(distVec)))/2;

%-------------------------------------------------------------------------------
% Convert to absolute correlations:
%-------------------------------------------------------------------------------
distVec0 = distVec; % keep distVec0 as the input distance matrix in all cases
% distVec changes only for abscorr options
if ismember(whatDistance,{'abscorr','abscorr_ii'})
    % Compute distances on absolute correlation distances (where sign of correlation
    % is irrelevant, it's the magnitude that's important):
    distVec = 1-abs(1-distVec);
end

% We need the matrix form also:
distMat = squareform(distVec);

%-------------------------------------------------------------------------------
% Do the linkage clustering:
%-------------------------------------------------------------------------------
fprintf(1,'Computing linkage information for %ux%u data using %s clustering...',...
            numItems,numItems,linkageMeth);
links = linkage(distVec,linkageMeth);
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% Compute the dendrogram
% ------------------------------------------------------------------------------
f = figure('color','w');
f.Position(3:4) = [1200,800];

% Get the dendrogram reordering:
ax1 = subplot(1,6,6);
ord = BF_linkageOrdering(distVec,links);
h_dend = dendrogram(links,0,'Orientation','right','Reorder',ord);
ax1.YDir = 'reverse'; % needs to be reversed to match the reversed y-axis of the imagesc plot

% Save the clustered distance matrix for output
distMat_cl = distMat(ord,ord);

% Make the dendrogram look aight
set(h_dend,'color','k','LineWidth',1)
ax1.YTickLabel = {};
xlabel('Distance');

%-------------------------------------------------------------------------------
% Cluster into groups:
%-------------------------------------------------------------------------------
% Cluster the dendrogram:
T = cluster(links,'cutoff',clusterThreshold,'criterion','distance');
numClusters = max(T);

fprintf(1,'Distance-based clustering with %u clusters\n',numClusters);

% Reorder members of each cluster by distance to other members of the cluster:
cluster_Groupi = cell(numClusters,1);
for j = 1:numClusters
    cluster_Groupi{j} = find(T==j);
    if length(T==j) > 1
        % Sort by increasing sum of distances to other members of the cluster
        [~,ix] = sort(sum(distMat(cluster_Groupi{j},cluster_Groupi{j})),'ascend');
        cluster_Groupi{j} = cluster_Groupi{j}(ix);
    end
end

% Reorder by decreasing cluster size
[~,ix] = sort(cellfun(@length,cluster_Groupi),'descend');
cluster_Groupi = cluster_Groupi(ix);

% Select the closest to cluster centre in each group
clusterCenters = cellfun(@(x)x(1),cluster_Groupi);

cluster_Groupi_cl = cellfun(@(x)arrayfun(@(y)find(ord==y),x),cluster_Groupi,'UniformOutput',0);
clusterCenters_cl = arrayfun(@(y)find(ord==y),clusterCenters);

% ------------------------------------------------------------------------------
% Now plot it:
% ------------------------------------------------------------------------------

% Make a square distance matrix, and plot as a pairwise similarity matrix
% (reordered by ord determined above)

ax2 = subplot(1,6,2:5);
switch whatDistance
case {'corr','abscorr_ii'}
    % Input is a correlation distance matrix, plot correlation coefficients:
    distMat0 = squareform(distVec0);
    BF_PlotCorrMat(1-distMat0(ord,ord),'auto');
case 'abscorr'
    distMat0 = squareform(distVec0);
    % Input is correlation distance matrix, plot absolute correlation coefficients:
    BF_PlotCorrMat(abs(1-distMat0(ord,ord)),'auto');
case 'general'
    % Input is a general distance matrix:
    BF_PlotCorrMat(distMat_cl);
otherwise
    warning('No special plotting options for distance metric: ''%s''',whatDistance);
    BF_PlotCorrMat(distMat_cl);
end

%-------------------------------------------------------------------------------
% Add rectangles to group highly correlated clusters:
%-------------------------------------------------------------------------------
rectangleColors = BF_getcmap('accent',5,1);
for i = 1:numClusters
    % Label cluster borders:
    rectangle('Position',[min(cluster_Groupi_cl{i})-0.5,min(cluster_Groupi_cl{i})-0.5, ...
                            length(cluster_Groupi_cl{i}),length(cluster_Groupi_cl{i})], ...
                    'EdgeColor',rectangleColors{1},'LineWidth',3);

    % Label cluster centers:
    rectangle('Position',[clusterCenters_cl(i)-0.5,clusterCenters_cl(i)-0.5,1,1], ...
                                'EdgeColor',rectangleColors{5},'LineWidth',3);
end

%-------------------------------------------------------------------------------
% Tweak plot details:
%-------------------------------------------------------------------------------

% Add object labels on y-axis if provided
if ~isempty(objectLabels)
    ax2.YTick = 1:numItems;
    ax2.YTickLabel = objectLabels(ord);
    ax2.TickLabelInterpreter = 'none';
else
    ax2.YTick = [];
end

% Remove X-ticks
ax2.XTick = [];

% Add a color bar:
cB = colorbar('southoutside');
switch whatDistance
case {'corr','abscorr_ii'}
    cB.Label.String = 'correlation coefficient';
case 'abscorr'
    cB.Label.String = 'absolute correlation coefficient';
case 'general'
    cB.Label.String = 'distance';
end

% Link axes:
linkaxes([ax1,ax2],'y');

% Fix limits
ax2.YLim = [0.5,numItems+0.5];

% Align with dendrogram
ax1.Position(1) = ax2.Position(1) + ax2.Position(3);
ax1.Position(4) = ax2.Position(4);
ax1.Position(2) = ax2.Position(2);

% Handles
if nargout > 3
    handles = struct();
    handles.f = f;
    handles.ax1 = ax1;
    handles.ax2 = ax2;
end

end
