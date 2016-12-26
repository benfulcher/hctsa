function TS_SimSearch(varargin)
% TS_SimSearch  Nearest neighbors of a given time series from an hctsa analysis.
%
% Nearest neighbors can provide a local context for a particular time series or
% operation.
%
%---INPUTS:
%
% targetID, the ID of the target time series or operation
% numNeighbors, the number of nearest neighbors to analyze
% whatPlots, a cell of plot types to show, e.g., {'matrix','network'}
%               (*) 'matrix' plots a pairwise similarity matrix
%               (*) 'scatter' plots scatter plots of outputs between target and neighbors
%               (*) 'network', a network visualization of neighbors
%
%---EXAMPLE USAGE:
%
% Find neighbors of time series (ID = 30), and visualize as a similarity matrix
% and network plot:
% TS_SimSearch(30,'whatPlots',{'matrix','network'})

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
% Check inputs, and set defaults
% ------------------------------------------------------------------------------

% Check inputs using the inputParser:
inputP = inputParser;

% targetID
default_targetID = 1;
check_targetID = @(x) validateattributes(x,{'numeric'},{'positive'});
addOptional(inputP,'targetID',default_targetID,check_targetID);

% tsOrOps
default_tsOrOps = 'ts';
valid_tsOrOps = {'ts','ops'};
check_tsOrOps = @(x) any(validatestring(x,valid_tsOrOps));
addOptional(inputP,'tsOrOps',default_tsOrOps,check_tsOrOps);

% whatDataFile
default_whatDataFile = 'norm';
check_whatDataFile = @(x)1;
addOptional(inputP,'whatDataFile',default_whatDataFile,check_whatDataFile);

% numNeighbors
default_numNeighbors = 20;
addOptional(inputP,'numNeighbors',default_numNeighbors,@isnumeric);

% whatPlots
default_whatPlots = {'matrix'};
check_whatPlots = @(x) iscell(x) || ischar(x);
addOptional(inputP,'whatPlots',default_whatPlots,check_whatPlots);

%-------------------------------------------------------------------------------
%% Parse inputs:
%-------------------------------------------------------------------------------
parse(inputP,varargin{:});

% Additional checks:
if strcmp(inputP.Results.whatDataFile,'cl')
    error('Using clustered permutations of data not suppored');
end

% Make variables from results of input parser:
targetID = inputP.Results.targetID;
tsOrOps = inputP.Results.tsOrOps;
numNeighbors = inputP.Results.numNeighbors;
whatDataFile = inputP.Results.whatDataFile;
whatPlots = inputP.Results.whatPlots;
clear inputP;

% ------------------------------------------------------------------------------
% Load data
% ------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatDataFile);
if strcmp(tsOrOps,'ts')
    dataStruct = TimeSeries;
    clear Operations
    clustStruct = TS_GetFromData(whatDataFile,'ts_clust');
    if isempty(clustStruct)
        % This should be set on normalization -- if missing for some reason, set as default now:
        clustStruct = struct('distanceMetric','none','Dij',[],...
                    'ord',1:size(TS_DataMat,1),'linkageMethod','none');
    end
else
    dataStruct = Operations;
    clear TimeSeries
    clustStruct = TS_GetFromData(whatDataFile,'op_clust');
    if isempty(clustStruct)
        % This should be set on normalization -- if missing for some reason, set as default now:
        clustStruct = struct('distanceMetric','none','Dij',[],...
                        'ord',1:size(TS_DataMat,2),'linkageMethod','none');
    end
end
numItems = length(dataStruct);

if numNeighbors == 0 % code for 'all'
    numNeighbors = numItems - 1;
else % specify a number of neighbours
    numNeighbors = min(numItems-1,numNeighbors);
end

% ------------------------------------------------------------------------------
% Match the specified index to the data structure
% ------------------------------------------------------------------------------
targetInd = find([dataStruct.ID]==targetID);
if isempty(targetInd)
    error('ID %u not found in the index for %s in %s.',targetID,tsOrOps,which(whatDataFile));
end

% ------------------------------------------------------------------------------
% Compute pairwise distances to all other objects
% ------------------------------------------------------------------------------
% (There is potential to store pairwise distance information in the HCTSA*.mat
% file for retrieval later). Use this if it exists, otherwise calculate for this one.

if isfield(clustStruct,'Dij') && ~isempty(clustStruct.Dij)
    % pairwise distances already computed, stored in the HCTSA .mat file
    fprintf(1,'Loaded %s distances from %s\n',clustStruct.distanceMetric,whatDataFile);
    Dij = squareform(clustStruct.Dij);
    Dj = Dij(:,targetInd);
else
    fprintf(1,'Computing %s distances to %u objects...',numItems);
    % Default: compute euclidean distances for time series and compute
    %           abs correlation distances for features
    switch tsOrOps
    case 'ts'
        Dj = bsxfun(@minus,TS_DataMat,TS_DataMat(targetInd,:));
        Dj = sqrt(mean(Dj.^2,2));
    case 'ops'
        % Is there a nicer way of computing abs correlations?
        Dj = zeros(numItems,1);
        for j = 1:numItems
            Dj(j) = 1 - abs(corr(TS_DataMat(:,targetInd),TS_DataMat(:,j)));
        end
    end
    fprintf(1,' Done.\n');
end

if strcmp(tsOrOps,'ops')
    % Now that distances are computed, transpose for consistency with later parts of the code
    TS_DataMat = TS_DataMat';
end

% ------------------------------------------------------------------------------
% Find N neighbors under the distance metric (used for Dj)
% ------------------------------------------------------------------------------

% Sort distances (ascending):
[~,dix] = sort(Dj,'ascend');

% Indices of nearest neighbors:
neighborInd = dix(1:numNeighbors+1);

% ------------------------------------------------------------------------------
% List matches to screen
% ------------------------------------------------------------------------------
fprintf(1,'\n---TARGET: %s---\n',dataStruct(targetInd).Name);
for j = 1:numNeighbors
    fprintf(1,'%u. %s (d = %.2f)\n',j,dataStruct(neighborInd(j+1)).Name,Dj(neighborInd(j+1)));
end
fprintf(1,'\n');

% ------------------------------------------------------------------------------
% Compute/retrieve pairwise distances between all neighbors
% ------------------------------------------------------------------------------
if isfield(clustStruct,'Dij') && ~isempty(clustStruct.Dij)
    % Use pre-computed distances:
    Dij = Dij(neighborInd,neighborInd)/sqrt(size(TS_DataMat,2)+1);
else
    % Recompute distances:
    Dij = squareform(pdist(TS_DataMat(neighborInd,:),'euclidean')/sqrt(size(TS_DataMat,2)+1));
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Plotting
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Scatter plot for top (up to 12)
% ------------------------------------------------------------------------------
if any(ismember(whatPlots,'scatter'))
    f = figure('color','w');
    for j = 1:min(12,numNeighbors)
        subplot(3,4,j);
        theNeighborInd = neighborInd(j+1); % Since exclude self-match (neighbor 1)
        plot(TS_DataMat(targetInd,:),TS_DataMat(theNeighborInd,:),'.k','MarkerSize',4)
        xlabel(sprintf('[%u] %s',dataStruct(targetInd).ID,dataStruct(targetInd).Name),'interpreter','none','FontSize',9);
        ylabel(sprintf('[%u] %s',dataStruct(theNeighborInd).ID,dataStruct(theNeighborInd).Name),'interpreter','none','FontSize',9);
        title(sprintf('Match %u. d = %.3f',j,Dj(theNeighborInd)));
        axis square
    end
    % Set width and height to make a reasonable size:
    f.Position = [f.Position(1:2),819,622];
end

% ------------------------------------------------------------------------------
% Clustered distance matrix
% ------------------------------------------------------------------------------
if any(ismember('matrix',whatPlots))
    % Do a quick cluster of distance matrix using euclidean distances, for visualization:
    ord = BF_ClusterReorder(Dij,'euclidean','average');
    Dij_clust = Dij(ord,ord);
    dataStruct_clust = dataStruct(neighborInd(ord));

    labels = cell(numNeighbors+1);
    for i = 1:numNeighbors+1
        labels{i} = sprintf('[%u] %s',dataStruct_clust(i).ID,dataStruct_clust(i).Name);
    end

    f = figure('color','w');

    % (I) Time-series annotations
    sp1 = subplot(1,5,1); ax1 = gca; box('on'); hold on
    ax1.YTick = (1:numNeighbors+1);
    ax1.YTickLabel = labels;
    ax1.YLim = [0.5,numNeighbors+1.5];
    tsLength = 100;
    ax1.XLim = [1,tsLength];
    xlabel('Time (samples)');
    ax1.TickLabelInterpreter = 'none';
    for j = 1:numNeighbors+1
        tsData = dataStruct(neighborInd(ord(j))).Data(1:tsLength);
        lengthHere = min(tsLength,length(tsData));
        plot(1:lengthHere,j-0.5+NormMinMax(tsData),'-k');
        if j < numNeighbors+1
            plot([1,tsLength],(j+0.5)*ones(2,1),':k')
        end
    end

    % (II) Pairwise similarity matrix
    sp2 = subplot(1,5,2:5); box('on'); ax2 = gca; hold on
    Dij_clust(logical(eye(numNeighbors+1))) = NaN; % zero diagonals mess things up
    Dij_clust(Dij_clust==0) = NaN; % all zeros mess things up
    dLims = [min(Dij_clust(~isnan(Dij_clust))),max(Dij_clust(~isnan(Dij_clust)))];
    imagesc(Dij_clust)
    if isfield(dataStruct,'Group')
        numGroups = length(unique([dataStruct_clust.Group]));
        dRescale = @(x) dLims(1) + numGroups/8*diff(dLims)*(-1 + 0.9999*(x - min(x))./(max(x)-min(x)));
        imagesc(0,1,dRescale([dataStruct_clust.Group]'))
        plot(ones(2,1)*0.5,[0.5,numNeighbors+1.5],'k')
        colormap([BF_getcmap('dark2',numGroups,0);BF_getcmap('redyellowblue',8,0)]);
        caxis([dLims(1)-diff(dLims)*numGroups/8,dLims(2)])
    else
        colormap(BF_getcmap('redyellowblue',8,0));
        caxis([min(Dij(Dij>0)),max(Dij(:))])
    end

    % Black rectangles over NaNs:
    [theNaNs_i,theNaNs_j] = find(isnan(Dij_clust));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end

    % Box the target:
    indexCl = find([dataStruct_clust.ID]==targetID);
    rectangle('Position',[indexCl-0.5,indexCl-0.5,1,1],'EdgeColor','w','LineWidth',2)
    plot(indexCl,indexCl,'*w')

    % Add a color bar:
    cB = colorbar('northoutside');
    cB.Label.String = 'Distance';
    if isfield(dataStruct,'Group')
        cB.Limits = dLims;
    end

    % Set axes:
    axis('square')
    if isfield(dataStruct,'Group')
        ax2.XLim = [-0.5,numNeighbors+1.5];
    else
        ax2.XLim = [0.5,numNeighbors+1.5];
    end
    ax2.XTick = 1:numNeighbors+1;
    ax2.XTickLabel = [dataStruct_clust.ID];

    ax2.YLim = [0.5,numNeighbors+1.5];
    ax2.YTick = 1:numNeighbors+1;
    ax2.YTickLabel = {};
    xlabel('ID');

    % Link axes:
    linkaxes([ax1,ax2],'y');

    % Reposition tight
    sp2.Position = [0.4,0.1,0.6,0.75];
    sp1.Position = [0.1,0.1,0.15,0.8];
    % Tight on left:
    sp1.Position(1) = sp2.Position(1) - sp1.Position(3) - 0.01;
    % Both have the same y and height:
    sp1.Position(2) = sp2.Position(2);
    sp1.Position(4) = sp2.Position(4);
end

% ------------------------------------------------------------------------------
% Network Visualization
% ------------------------------------------------------------------------------
if any(ismember(whatPlots,'network'))

    % First ensure that no more than a maximum number of neighbors are plotted:
    numNetwork = 35; % max neighbors for network vis
    if numNeighbors > numNetwork
        warning('Displaying a maximum of %u neighbors for the network visualization',numNetwork)
    end
    numNetwork = min(numNeighbors,numNetwork); % number neighbours for network

    % Define the adjacency matrix, A (target node is 1)
    A = 1 - Dij(1:numNetwork,1:numNetwork);

    % Remove the diagonal:
    A(logical(eye(size(A)))) = 0;

    % Assign group labels (or just distinguish the target)
    if isfield(dataStruct,'Group')
        nodeLabels = [dataStruct.Group];
        nodeLabels = nodeLabels(dix(1:numNetwork));
        nodeLabels(1) = max(nodeLabels) + 1;
    else
        nodeLabels = 2*ones(numNetwork,1);
        nodeLabels(1) = 1;
    end

    if strcmp(tsOrOps,'ts')
        dataLabels = {dataStruct(dix(1:numNetwork)).Data};
    else
        dataLabels = {};
    end

    NetVis_netvis(A,'k',0.01,'textLabels',{dataStruct(dix(1:numNetwork)).Name},...
                    'linkThresh',[0.9,0.8,0.7,0.6],...
                    'nodeLabels',nodeLabels,...
                    'dataLabels',dataLabels,...
                    'colorMap','set1');
end

% ------------------------------------------------------------------------------
function X = NormMinMax(x)
    X = (x-min(x))/(max(x)-min(x));
end
% ------------------------------------------------------------------------------

end
