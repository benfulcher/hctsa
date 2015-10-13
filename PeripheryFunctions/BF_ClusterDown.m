function [distMat_cl,cluster_Groupi,ord] = BF_ClusterDown(distMat,numClusters,varargin)
% BF_ClusterDown    Reduce a pairwise similarity matrix into smaller clusters.
%
% Yields a visualization of a pairwise similarity matrix, with an attempt to
% deduce a set of smaller clusters of objects showing highly correlated behavior.

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

%-------------------------------------------------------------------------------
% Check Inputs
%-------------------------------------------------------------------------------
% numClusters -- the number of clusters to generate from your data:
if nargin < 2 || isempty(numClusters)
    numClusters = 5; % removes very highly-correlated operations
end

% Custom plotting options using an input parser
inputP = inputParser;

% Object labels for axis labels
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
objectLabels = inputP.Results.objectLabels;
whatDistance = inputP.Results.whatDistance;
% errTh = inputP.Results.errTh;
% miOrCorr = inputP.Results.miOrCorr;
% clusterMeth = inputP.Results.clusterMeth;
linkageMeth = inputP.Results.linkageMeth;
% plotBig = inputP.Results.plotBig;
clear inputP;

%-------------------------------------------------------------------------------
% Check squareform:
%-------------------------------------------------------------------------------
if any(size(distMat)==1)
    % A vector, probably a squareform version
    distMat = squareform(distMat);
end

%-------------------------------------------------------------------------------
% Do the linkage clustering:
%-------------------------------------------------------------------------------
numItems = length(distMat);
fprintf(1,'Computing linkage information for %ux%u data using %s clustering...',...
            numItems,numItems,linkageMeth);
links = linkage(distMat,linkageMeth);
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% Compute the dendrogram
% ------------------------------------------------------------------------------

f = figure('color','w');
f.Position(3) = 1200;
f.Position(4) = 800;

% Get the dendrogram reordering:
subplot(1,6,6); ax1 = gca;
[h_dend,~,ord] = dendrogram(links,0,'Orientation','right');
ord = fliplr(ord); % so everything matches up with the dendrogram
% Reorder the distance matrix by dendrogram ordering: [could add optimalleaforder]
distMat_cl = distMat(ord,ord);

% Make the dendrogram look aight
set(h_dend,'color','k','LineWidth',1)
ax1.YTickLabel = {};
xlabel('Distance');

% Compute clustering into groups for many different cluster numbers:
numClusterings = length(numClusters);
cluster_Groupi = cell(numClusterings,1);
cluster_Groupi_cl = cell(numClusterings,1);
clusterCenters_cl = cell(numClusterings,1);
clusterCenters = cell(numClusterings,1);

for i = 1:numClusterings
    fprintf(1,'Distance-based clustering with %u clusters\n',numClusters(i))

    % Cluster the dendrogram:
    T = cluster(links,'maxclust',numClusters(i));
    numClusters(i) = max(T);

    % Reorder members of each cluster by distance to other members of the cluster:
    cluster_Groupi{i} = cell(numClusters(i),1);
    for j = 1:numClusters(i)
        cluster_Groupi{i}{j} = find(T==j);
        if length(T==j) > 1
            % Sort by increasing sum of distances to other members of the cluster
            [~,ix] = sort(sum(distMat(cluster_Groupi{i}{j},cluster_Groupi{i}{j})),'ascend');
            cluster_Groupi{i}{j} = cluster_Groupi{i}{j}(ix);
        end
    end

    % Reorder by decreasing cluster size
    [~,ix] = sort(cellfun(@length,cluster_Groupi{i}),'descend');
    cluster_Groupi{i} = cluster_Groupi{i}(ix);

    % Select the closest to cluster centre in each group
    clusterCenters{i} = cellfun(@(x)x(1),cluster_Groupi{i});

    cluster_Groupi_cl{i} = cellfun(@(x)arrayfun(@(y)find(ord==y),x),cluster_Groupi{i},'UniformOutput',0);
    clusterCenters_cl{i} = arrayfun(@(y)find(ord==y),clusterCenters{i});
end


% ------------------------------------------------------------------------------
% Now plot it:
% ------------------------------------------------------------------------------
% Pick a given clustering and plot it
cLevel = min(1,numClusterings); % plot the first clustering

% Plot as a similarity matrix:
subplot(1,6,2:5)
ax2 = gca;
switch whatDistance
case 'corr'
    % Input is a correlation matrix:
    BF_PlotCorrMat(1-distMat_cl);
case 'abscorr'
    % Input is a correlation matrix -> absolute value distances
    BF_PlotCorrMat(1-abs(distMat_cl));
case 'general'
    % Input is a general distance matrix:
    BF_PlotCorrMat(distMat_cl);
otherwise
    error('Unknown distance ''%s''',whatDistance);
end

% Add rectangles to indicate highly correlated clusters of statistics:
rectangleColors = BF_getcmap('accent',5,1);
for i = 1:numClusters
    % Cluster borders:
    rectangle('Position',[min(cluster_Groupi_cl{cLevel}{i})-0.5,min(cluster_Groupi_cl{cLevel}{i})-0.5, ...
                            length(cluster_Groupi_cl{cLevel}{i}),length(cluster_Groupi_cl{cLevel}{i})], ...
                    'EdgeColor',rectangleColors{1},'LineWidth',3);

    % Cluster centers:
    rectangle('Position',[clusterCenters_cl{cLevel}(i)-0.5,clusterCenters_cl{cLevel}(i)-0.5,1,1], ...
                                'EdgeColor',rectangleColors{5},'LineWidth',3);
end

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
case 'corr'
    cB.Label.String = 'correlation coefficient';
case 'general'
    cB.Label.String = 'distance';
end

% Link axes:
linkaxes([ax1,ax2],'y');

% Fix limits
ax2.YLim = [0.5,numItems+0.5];

% Align with dendrogram
ax1.Position(1) = ax2.Position(1) + ax2.Position(3);
ax1.Position(2) = ax2.Position(2);
ax1.Position(4) = ax2.Position(4);

end
