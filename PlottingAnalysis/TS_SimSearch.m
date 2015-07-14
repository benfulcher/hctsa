% ------------------------------------------------------------------------------
% TS_SimSearch
% ------------------------------------------------------------------------------
%
% Finds nearest neighbors of a given item, providing a local context for a
% particular time series or operation.
%
%---INPUTS:
%
%
% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TS_SimSearch(tsOrOps,whatDataFile,targetID,numNeighbors)

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(tsOrOps)
    tsOrOps = 'ts';
    fprintf(1,'Time series by default.\n');
end
if ~ismember(tsOrOps,{'ts','ops'})
    error('Unknown specifier ''%s''. Should be ''ts'' or ''ops''.',tsOrOps);
end

if nargin < 2
    whatDataFile = 'norm';
end

if nargin < 3
    targetID = 1;
end

if nargin < 4 || isempty(numNeighbors)
    numNeighbors = 20;
end

% ------------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------------

[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatDataFile);
if strcmp(tsOrOps,'ts')
    dataStruct = TimeSeries;
    clear Operations
    tmp = load(whatDataFile,'ts_clust');
    clustStruct = tmp.ts_clust; clear tmp
else
    dataStruct = Operations;
    clear TimeSeries
    tmp = load(whatDataFile,'op_clust');
    clustStruct = tmp.op_clust; clear tmp
end
numItems = length(dataStruct);
numNeighbors = min(numItems,numNeighbors) - 1;

% ------------------------------------------------------------------------------
% Match the specified index to the data structure
% ------------------------------------------------------------------------------
targetInd = find([dataStruct.ID]==targetID);
if isempty(targetInd)
    error('ID %u not found in the index for %s in %s.',targetID,dataStruct,which(theFile));
end

% ------------------------------------------------------------------------------
% Compute pairwise distances to all other objects
% ------------------------------------------------------------------------------
% (There is potential to store pairwise distance information in the HCTSA_loc
% file for retrieval later). I guess use this if it exists, otherwise just
% calculate for this one.


if isfield(clustStruct,'Dij')
    % pairwise distances already computed, stored in the HCTSA .mat file
    fprintf(1,'Loaded %s distances from %s\n',clustStruct.distanceMetric,whatDataFile)
    Dij = squareform(clustStruct.Dij);
    Dj = Dij(:,targetInd);
    clear Dij;
else
    fprintf(1,'Computing distances to %u objects...',numItems);
    switch tsOrOps
    case 'ts'
        Dj = bsxfun(@minus,TS_DataMat,TS_DataMat(targetInd,:));
        Dj = sqrt(mean(Dj.^2,2));
    case 'ops'
        % Is there a nicer way of computing correlations?
        Dj = zeros(numItems,1);
        for j = 1:numItems
            Dj(j) = corr(TS_DataMat(:,targetInd),TS_DataMat(:,j));
        end
    end
    fprintf(1,' Done.\n');
end

if strcmp(tsOrOps,'ops')
    % Now that distances are computed, transpose for consistency with later parts of the code
    TS_DataMat = TS_DataMat';
end

% ------------------------------------------------------------------------------
% Find N neighbors under a given distance metric
% ------------------------------------------------------------------------------

% Sort distances (ascending):
[~,dix] = sort(Dj,'ascend');

% Indices of nearest neighbors:
neighborInd = dix(1:numNeighbors+1);

% ------------------------------------------------------------------------------
% List matches to screen
% ------------------------------------------------------------------------------
fprintf(1,'---%s---\n',dataStruct(targetInd).Name);
for j = 2:numNeighbors+1
    fprintf(1,'%u. %s (d = %.2f)\n',j-1,dataStruct(neighborInd(j)).Name,Dj(neighborInd(j)));
end

% ------------------------------------------------------------------------------
% Compute pairwise distances between all neighbors
% ------------------------------------------------------------------------------
Dij = squareform(pdist(TS_DataMat(neighborInd,:),'euclidean')/sqrt(size(TS_DataMat,2)+1));

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Plotting
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Scatter plot for top (up to 12)
% ------------------------------------------------------------------------------
f = figure('color','w');
for j = 2:min(12,numNeighbors)+1
    subplot(3,4,j-1);
    plot(TS_DataMat(targetInd,:),TS_DataMat(neighborInd(j),:),'.k','MarkerSize',4)
    xlabel(sprintf('[%u] %s',dataStruct(targetInd).ID,dataStruct(targetInd).Name),'interpreter','none','FontSize',9);
    ylabel(sprintf('[%u] %s',dataStruct(neighborInd(j)).ID,dataStruct(neighborInd(j)).Name),'interpreter','none','FontSize',9);
    title(sprintf('Match %u. d = %.3f',j,Dj(neighborInd(j))))
    axis square
end

% ------------------------------------------------------------------------------
% Clustered distance matrix
% ------------------------------------------------------------------------------
% Do a quick cluster of distance matrix using euclidean distances, for visualization:
ord = BF_ClusterReorder(Dij,'euclidean','average');
Dij_clust = Dij(ord,ord);
dataStruct_clust = dataStruct(neighborInd(ord));

labels = cell(numNeighbors+1);
for i = 1:numNeighbors+1
    labels{i} = sprintf('[%u] %s',dataStruct_clust(i).ID,dataStruct_clust(i).Name);
end

f = figure('color','w');
sp1 = subplot(1,5,1); ax = gca; hold on
ax.YTick = (1:numNeighbors+1);
ax.YTickLabel = labels;
ax.YLim = [0.5,numNeighbors+1.5];
tsLength = 100;
ax.XLim = [1,tsLength];
xlabel('Time (samples)');
ax.TickLabelInterpreter = 'none';
for j = 1:numNeighbors+1
    tsData = dataStruct(neighborInd(ord(j))).Data(1:100);
    lengthHere = min(tsLength,length(tsData));
    plot(1:lengthHere,j-0.5+NormMinMax(tsData),'-k');
end
sp2 = subplot(1,5,2:5); ax = gca; hold on
imagesc(Dij_clust)
caxis([min(Dij(Dij>0)),max(Dij(:))])
colormap(BF_getcmap('redyellowblue',8,0));

% Add a color bar:
cB = colorbar('northoutside');
cB.Label.String = 'Distance';

% Box the target:
indexCl = find([dataStruct_clust.ID]==targetID);
rectangle('Position',[indexCl-0.5,indexCl-0.5,1,1])
plot(indexCl,indexCl,'*k')

% Set axes:
axis('square')
ax.XLim = [0.5,numNeighbors+1.5];
ax.YLim = [0.5,numNeighbors+1.5];
ax.XTick = 1:numNeighbors+1;
ax.XTickLabel = [dataStruct_clust.ID];
ax.YTick = 1:numNeighbors+1;
ax.YTickLabel = {};
xlabel('ID');

% Reposition tight
sp2.Position = [0.4,0.1,0.6,0.75];
sp1.Position = [0.1,0.1,0.15,0.8];
% Tight on left:
sp1.Position(1) = sp2.Position(1) - sp1.Position(3) - 0.01;
% Both have the same y and height:
sp1.Position(2) = sp2.Position(2);
sp1.Position(4) = sp2.Position(4);

% ------------------------------------------------------------------------------
% Network Visualization
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Make a network of neighbours
% ------------------------------------------------------------------------------
numNetwork = 30; % max neighbors for network vis
numNetwork = min(numNeighbors,numNetwork); % number neighbours for network

% Acl of distances, and a labelscl of labels
A = 1 - Dij(1:numNetwork,1:numNetwork);

if isfield(dataStruct,'Group')
    nodeLabels = [dataStruct.Group];
    nodeLabels = nodeLabels(dix(1:numNetwork));
    nodeLabels(1) = max(nodeLabels) + 1;
else
    nodeLabels = 2*ones(numNetwork,1);
    nodeLabels(1) = 1;
end

A(logical(eye(size(A)))) = 0;
NetVis_netvis(A,0.01,{dataStruct.Name},[0.9,0.8,0.7,0.6],nodeLabels,[],[1000,1],'set1');

% ------------------------------------------------------------------------------
function X = NormMinMax(x)
    X = (x-min(x))/(max(x)-min(x));
end
% ------------------------------------------------------------------------------

end
