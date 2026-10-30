function TS_SimSearch(varargin)
% TS_SimSearch  Nearest neighbors of a given time series from an hctsa analysis.
%
% Nearest neighbors can provide a local context for a particular time series or
% operation.
%
%---(OPTIONAL) INPUTS:
% targetID, the ID of the target time series or operation
% tsOrOps, 'ts' (for time series) or 'ops' for operations (features)
% numNeighbors, the number of nearest neighbors to analyze
% whatPlots, a cell of plot types to show, e.g., {'matrix','network'}
%               (*) 'matrix' plots a pairwise similarity matrix
%               (*) 'scatter' plots scatter plots of outputs between target and neighbors
%               (*) 'network', a network visualization of neighbors
% whatDistMetric, the distance metric used to find neighbors, when distances are
%               not stored in the data (cf. TS_Cluster): 'Euclidean', 'Spearman'
%               (1 - Spearman correlation), 'Pearson' or 'corr' (1 - Pearson
%               correlation), 'abscorr' (1 - |Spearman correlation|), or 'corr_fast'.
%               Defaults match TS_Cluster: 'Euclidean' for time series and
%               'corr_fast' for operations. When distances are stored in the data,
%               those are used (and the metric they were computed with).
%               Distances are defined as in BF_pdist, so the neighbors are the same
%               whether or not TS_Cluster has been run.
%
%---EXAMPLE USAGE:
% Find neighbors of time series (ID = 30), and visualize as a similarity matrix
% and network plot:
% TS_SimSearch('targetID',30,'whatPlots',{'matrix','network'})

%-------------------------------------------------------------------------------
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
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs, and set defaults
%-------------------------------------------------------------------------------

% Check inputs using the inputParser:
inputP = inputParser;

% targetID
default_targetID = 1;
check_targetID = @(x) validateattributes(x,{'numeric'},{'positive'});
addParameter(inputP,'targetID',default_targetID,check_targetID);

% tsOrOps
default_tsOrOps = 'ts';
valid_tsOrOps = {'ts','ops'};
check_tsOrOps = @(x) any(validatestring(x,valid_tsOrOps));
addParameter(inputP,'tsOrOps',default_tsOrOps,check_tsOrOps);

% whatData
default_whatData = 'norm';
check_whatData = @(x)true;
addParameter(inputP,'whatData',default_whatData,check_whatData);

% numNeighbors
default_numNeighbors = 20;
addParameter(inputP,'numNeighbors',default_numNeighbors,@isnumeric);

% whatPlots
default_whatPlots = {'matrix'};
check_whatPlots = @(x) iscell(x) || ischar(x);
addParameter(inputP,'whatPlots',default_whatPlots,check_whatPlots);

% whatDistMetric
default_whatDistMetric = '';
check_whatDistMetric = @(x) ischar(x) && (isempty(x) || ...
        ismember(x,{'Euclidean','Spearman','Pearson','corr','abscorr','corr_fast'}));
addParameter(inputP,'whatDistMetric',default_whatDistMetric,check_whatDistMetric);

%-------------------------------------------------------------------------------
%% Parse inputs:
%-------------------------------------------------------------------------------
parse(inputP,varargin{:});

% Additional checks:
if strcmp(inputP.Results.whatData,'cl')
    error('Using clustered permutations of data not supported');
end

% Make variables from results of input parser:
targetID = inputP.Results.targetID;
tsOrOps = inputP.Results.tsOrOps;
numNeighbors = inputP.Results.numNeighbors;
whatData = inputP.Results.whatData;
whatPlots = inputP.Results.whatPlots;
if ischar(whatPlots)
    whatPlots = {whatPlots};
end
whatDistMetric = inputP.Results.whatDistMetric;
clear inputP;

%-------------------------------------------------------------------------------
% Load data
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatData] = TS_LoadData(whatData);
switch tsOrOps
case 'ts'
    dataTable = TimeSeries;
    clustStruct = TS_GetFromData(whatData,'ts_clust');
case 'ops'
    dataTable = Operations;
    clustStruct = TS_GetFromData(whatData,'op_clust');
    % Transpose for consistency with later parts of the code
    % (items are rows)
    TS_DataMat = TS_DataMat';
end

if isempty(clustStruct)
    % This should be set on normalization -- if missing for some reason, set as default now:
    clustStruct = struct('distanceMetric','none','Dij',[],...
                'ord',1:size(TS_DataMat,1),'linkageMethod','none');
end
numItems = height(dataTable);
clear('TimeSeries','Operations')

% Distances: use those stored in the data file (by TS_Cluster) if they exist,
% otherwise compute them
hasStoredDist = isfield(clustStruct,'Dij') && ~isempty(clustStruct.Dij);
if hasStoredDist
    % The stored distances define the metric
    distMetric = clustStruct.distanceMetric;
    if ~isempty(whatDistMetric)
        fprintf(1,'Using the stored %s distances (ignoring whatDistMetric = ''%s'').\n',...
                    distMetric,whatDistMetric);
    end
else
    if isempty(whatDistMetric)
        switch tsOrOps
        case 'ts'
            whatDistMetric = 'Euclidean';
        case 'ops'
            whatDistMetric = 'corr_fast';
        end
        fprintf(1,'Using default distance metric: %s\n',whatDistMetric);
    end
    % The name of the metric in BF_pdist:
    switch whatDistMetric
    case 'Euclidean'
        distMetric = 'euclidean';
    case 'corr'
        distMetric = 'Pearson';
    otherwise
        distMetric = whatDistMetric;
    end
end
% Is it a correlation distance (1 - rho or 1 - |rho|)?
if ismember(distMetric,{'abscorr','abscorr_fast'})
    corrKind = 'abs';
elseif ismember(distMetric,{'Spearman','Pearson','PearsonCorr','corr','correlation','corr_fast'})
    corrKind = 'signed';
else
    corrKind = 'none';
end

if numNeighbors == 0 % code for 'all'
    numNeighbors = numItems - 1;
else % specify a number of neighbours
    numNeighbors = min(numItems-1,numNeighbors);
end

%-------------------------------------------------------------------------------
% Match the specified index to the data structure
%-------------------------------------------------------------------------------
targetInd = find(dataTable.ID==targetID);
if isempty(targetInd)
    error('ID %u not found in the index for %s in %s.',targetID,tsOrOps,which(whatData));
else
    fprintf(1,'\n---TARGET: [%u] %s (%s)---\n',...
        dataTable.ID(targetInd),dataTable.Name{targetInd},dataTable.Keywords{targetInd});
end

%-------------------------------------------------------------------------------
% Compute distance from target to all other objects
%-------------------------------------------------------------------------------
% (There is potential to store pairwise distance information in the HCTSA*.mat
% file for retrieval later). Use this if it exists, otherwise calculate for this one.

if hasStoredDist
    % pairwise distances already computed, stored in the HCTSA .mat file
    fprintf(1,'Loaded %s distances from %s\n',distMetric,whatData);
    Dij = squareform(clustStruct.Dij);
    Dj = Dij(:,targetInd);
else
    % Compute distances, as defined by BF_pdist:
    % (Note that TS_DataMat has been transposed in the case of 'ops')
    fprintf(1,'Computing %s distances to %u other objects...',distMetric,numItems-1);
    Dj = TargetDistances(TS_DataMat,targetInd,distMetric);
    fprintf(1,' Done.\n');
end

%-------------------------------------------------------------------------------
% Find N neighbors under the distance metric (used for Dj)
%-------------------------------------------------------------------------------
% Sort distances (ascending), putting the target first by index (other items can have
% distance zero to it, e.g., duplicates, and its own distance is undefined if it has NaNs):
Dsort = Dj;
Dsort(targetInd) = -Inf;
[~,dix] = sort(Dsort,'ascend');

% Indices of the target and its nearest neighbors:
neighborInd = dix(1:numNeighbors+1);

%-------------------------------------------------------------------------------
% List matches to screen
%-------------------------------------------------------------------------------
for j = 1:numNeighbors
    theInd = neighborInd(j+1);
    switch corrKind
    case 'signed'
        fprintf(1,'%u. [%u] %s (%s): rho = %.2f\n',j,dataTable.ID(theInd),...
                dataTable.Name{theInd},dataTable.Keywords{theInd},1 - Dj(theInd));
    case 'abs'
        fprintf(1,'%u. [%u] %s (%s): |rho| = %.2f\n',j,dataTable.ID(theInd),...
                dataTable.Name{theInd},dataTable.Keywords{theInd},1 - Dj(theInd));
    otherwise
        fprintf(1,'%u. [%u] %s (%s): d = %.2f\n',j,dataTable.ID(theInd),...
                dataTable.Name{theInd},dataTable.Keywords{theInd},Dj(theInd));
    end
end
fprintf(1,'\n');

%-------------------------------------------------------------------------------
% Compute/retrieve pairwise distances between all neighbors
%   (needed for 'matrix' and 'network')
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'matrix')) || any(ismember(whatPlots,'network'))
    if hasStoredDist
        % Use pre-computed distances:
        Dij = Dij(neighborInd,neighborInd);
    else
        % Compute distances (on the same scale as the stored ones):
        Dij = BF_pdist(TS_DataMat(neighborInd,:),distMetric,false,[],true);
    end
end

%===============================================================================
%===============================================================================
% Plotting
%===============================================================================
%===============================================================================

%-------------------------------------------------------------------------------
% Scatter plot for top (up to 12)
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'scatter'))
    doRank = false;
    f = figure('color','w');
    numScatters = min(12,numNeighbors);
    h = cell(numScatters,1);

    xData = TS_DataMat(targetInd,:);
    [~,ix] = sort(xData,'ascend');
    xDataRanks = 1:length(xData);
    xDataRanks(ix) = xDataRanks;

    for j = 1:numScatters
        h{j} = subplot(3,4,j);
        % Get data and rank it:
        theNeighborInd = neighborInd(j+1); % Since exclude self-match (neighbor 1)
        yData = TS_DataMat(theNeighborInd,:);
        [~,iy] = sort(yData,'ascend');
        yDataRanks = 1:length(yData);
        yDataRanks(iy) = yDataRanks;

        hold('on')
        h{j}.UserData = struct();
        h{j}.UserData.hRaw = plot(xData,yData,'.k','MarkerSize',4,'Visible','off');
        h{j}.UserData.hRank = plot(xDataRanks,yDataRanks,'.b','MarkerSize',4,'Visible','off');
        if doRank
            h{j}.UserData.hRank.Visible = 'on'; %#ok<UNRCH>
        else
            h{j}.UserData.hRaw.Visible = 'on';
        end
        xlabel(sprintf('[%u] %s',dataTable.ID(targetInd),dataTable.Name{targetInd}),...
                                        'interpreter','none','FontSize',9);
        ylabel(sprintf('[%u] %s',dataTable.ID(theNeighborInd),dataTable.Name{theNeighborInd}),...
                                        'interpreter','none','FontSize',9);
        title(sprintf('Match %u: d = %.3f',j,Dj(theNeighborInd)));
        axis('square')

        h{j}.ButtonDownFcn = {@swapRanked};
        h{j}.HitTest = 'on';
        h{j}.UserData.doRank = doRank;
        % h{j}.UserData.hRaw = hRaw;
        % h{j}.UserData.hRank = hRank;
        % h{j}.UserData.yData = yData;
        % h{j}.UserData.xDataRanks = xDataRanks;
        % h{j}.UserData.yDataRanks = yDataRanks;
    end
    % Set width and height to make a reasonable size:
    f.Position = [f.Position(1:2),819,622];

    fprintf(1,'DID YOU KNOW YOU CAN CLICK ON THE SCATTER PLOT TO SWITCH BETWEEN PLOTTING VALUES/RANKS?\n');
    fprintf(1,'GIVE IT A GO?!\n');
end

%-------------------------------------------------------------------------------
% Clustered distance matrix
%-------------------------------------------------------------------------------
if any(ismember('matrix',whatPlots))
    % Do a quick cluster of distance matrix using euclidean distances
    % (just for visualization):
    ord = BF_ClusterReorder(Dij,'euclidean','average');
    % BF_ClusterReorder leaves out items with NaN distances; show them at the end:
    ord = [ord(:);setdiff((1:numNeighbors+1)',ord(:))];
    Dij_clust = Dij(ord,ord);
    dataTable_clust = dataTable(neighborInd(ord),:);

    labels = cell(numNeighbors+1,1);
    for i = 1:numNeighbors+1
        labels{i} = sprintf('[%u] %s',dataTable_clust.ID(i),dataTable_clust.Name{i});
    end

    figure('color','w');

    switch tsOrOps
    case 'ts'
        % (I) Time-series annotations
        ax1 = subplot(1,5,1); box('on'); hold on
        ax1.YTick = (1:numNeighbors+1);
        ax1.YTickLabel = labels;
        ax1.YLim = [0.5,numNeighbors+1.5];
        tsLength = 100;
        ax1.XLim = [1,tsLength];
        xlabel('Time (samples)');
        ax1.TickLabelInterpreter = 'none';
        NormMinMax = @(x) (x-min(x))/(max(x)-min(x));
        for j = 1:numNeighbors+1
            tsData = dataTable.Data{neighborInd(ord(j))};
            lengthHere = min(tsLength,length(tsData));
            plot(1:lengthHere,j-0.5+NormMinMax(tsData(1:lengthHere)),'-k');
            if j < numNeighbors+1
                plot([1,tsLength],(j+0.5)*ones(2,1),':k')
            end
        end
        % (II) Pairwise similarity matrix
        ax2 = subplot(1,5,2:5);
        box('on'); hold('on')
    case 'ops'
        % Pairwise similarity matrix only
        ax2 = gca;
        box('on'); hold('on')
    end

    Dij_clust(logical(eye(numNeighbors+1))) = NaN; % zero diagonals mess things up
    Dij_clust(Dij_clust==0) = NaN; % all zeros mess things up
    dLims = [min(Dij_clust(~isnan(Dij_clust))),max(Dij_clust(~isnan(Dij_clust)))];
    if isempty(dLims) || dLims(2) <= dLims(1)
        % (no distinct nonzero distances to display, e.g., perfectly correlated neighbors)
        dLims = [0,1];
    end
    imagesc(Dij_clust)

    % Get group labeling if possible (time series):
    [~,~,groupLabelsInteger,numClasses] = TS_ExtractGroupLabels(dataTable);

    % Add class labels as rectangles along the edge…
    if numClasses > 1
        % The class of each displayed item (unlabeled items, coded 0, are left blank),
        % scaled so that classes 1,...,numClasses take the class colors of the colormap:
        groupLabelsShown = groupLabelsInteger(neighborInd(ord));
        groupLabelsShown(groupLabelsShown==0) = NaN;
        dRescale = @(x) dLims(1) + numClasses/8*diff(dLims)*(-1 + 0.9999*(x - 1)./(numClasses - 1));
        imagesc(0,1,dRescale(groupLabelsShown))
        plot(ones(2,1)*0.5,[0.5,numNeighbors+1.5],'k')
        colormap([BF_GetColorMap('dark2',numClasses,0);BF_GetColorMap('redyellowblue',8,0)]);
        clim([dLims(1)-diff(dLims)*numClasses/8,dLims(2)])
    else
        colormap(BF_GetColorMap('redyellowblue',8,0));
        clim(dLims)
    end

    % Black rectangles over NaNs:
    [theNaNs_i,theNaNs_j] = find(isnan(Dij_clust));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end

    % Box the target:
    indexCl = find(dataTable_clust.ID==targetID);
    rectangle('Position',[indexCl-0.5,indexCl-0.5,1,1],'EdgeColor','w','LineWidth',2)
    plot(indexCl,indexCl,'*w')

    % Add a color bar:
    cB = colorbar('northoutside');
    cB.Label.String = sprintf('%s distance',distMetric);
    if numClasses > 1
        cB.Limits = dLims;
    end

    % Set axes:
    axis('square')
    if numClasses > 1
        ax2.XLim = [-0.5,numNeighbors+1.5];
    else
        ax2.XLim = [0.5,numNeighbors+1.5];
    end
    ax2.XTick = 1:numNeighbors+1;
    ax2.XTickLabel = dataTable_clust.ID;
    ax2.XTickLabelRotation = 60;
    ax2.FontSize = 9;
    xlabel('ID');

    ax2.YLim = [0.5,numNeighbors+1.5];
    ax2.YTick = 1:numNeighbors+1;
    switch tsOrOps
    case 'ts'
        ax2.YTickLabel = [];
    case 'ops'
        ax2.YTickLabel = labels;
        ax2.TickLabelInterpreter = 'none';
        ylabel('Name');
    end

    % Link axes:
    if strcmp(tsOrOps,'ts')
        linkaxes([ax1,ax2],'y');
        % Reposition tight
        ax2.Position = [0.4,0.1,0.6,0.75];
        ax1.Position = [0.1,0.1,0.15,0.8];
        % Tight on left:
        ax1.Position(1) = ax2.Position(1) - ax1.Position(3) - 0.01;
        % Both have the same y and height:
        ax1.Position(2) = ax2.Position(2);
        ax1.Position(4) = ax2.Position(4);
    end
end

%-------------------------------------------------------------------------------
% Network Visualization
%-------------------------------------------------------------------------------
if any(ismember(whatPlots,'network'))

    % First ensure that no more than a maximum number of neighbors are plotted:
    numNetwork = 35; % max neighbors for network vis
    if numNeighbors > numNetwork
        warning('Displaying a maximum of %u neighbors for the network visualization',numNetwork)
    end
    numNetwork = min(numNeighbors,numNetwork) + 1; % number of nodes: the target and its neighbors

    % Define the adjacency matrix, A (target node is 1), as a similarity in [0,1]:
    Dnet = Dij(1:numNetwork,1:numNetwork);
    if strcmp(corrKind,'none')
        % Relative to the largest distance displayed
        maxDist = max(Dnet(:),[],'omitnan');
        if isempty(maxDist) || isnan(maxDist) || maxDist==0
            maxDist = 1;
        end
        A = 1 - Dnet/maxDist;
    else
        % Correlation similarity (negative correlations are not linked)
        A = max(1 - Dnet,0);
    end
    A(isnan(A)) = 0;

    % Remove the diagonal:
    A(logical(eye(size(A)))) = 0;

    % Assign group labels (or just distinguish the target)
    if ismember('Group',dataTable.Properties.VariableNames)
        % Integer group labels, with the target as its own group:
        [~,~,nodeLabels] = TS_ExtractGroupLabels(dataTable);
        nodeLabels = nodeLabels(neighborInd(1:numNetwork));
        nodeLabels(1) = max(nodeLabels) + 1;
        % Relabel as 1,2,...,numGroups for the network visualization:
        [~,~,nodeLabels] = unique(nodeLabels);
    else
        nodeLabels = 2*ones(numNetwork,1);
        nodeLabels(1) = 1;
    end

    if strcmp(tsOrOps,'ts')
        dataLabels = dataTable.Data(neighborInd(1:numNetwork));
    else
        dataLabels = {};
    end

    % (links are classed by quantiles of A, so the legend gives the strongest x% of entries)
    NetVis_netvis(A,'k',0.01,'textLabels',dataTable.Name(neighborInd(1:numNetwork)),...
                    'linkThresh',[0.9,0.8,0.7,0.6],...
                    'thresholdByProportion',true,...
                    'nodeLabels',nodeLabels,...
                    'dataLabels',dataLabels,...
                    'colorMap','set1',...
                    'makeFigure',true);
end


%-------------------------------------------------------------------------------
function swapRanked(hObject,~)
    % Swap between raw and ranked values for a scatter plot
    if hObject.UserData.doRank
        hObject.UserData.hRaw.Visible = 'on';
        hObject.UserData.hRank.Visible = 'off';
        hObject.XLabel.String = hObject.XLabel.String(1:end-8);
        hObject.YLabel.String = hObject.YLabel.String(1:end-8);
    else
        hObject.UserData.hRaw.Visible = 'off';
        hObject.UserData.hRank.Visible = 'on';
        hObject.XLabel.String = sprintf('%s%s',hObject.XLabel.String,' [RANKS]');
        hObject.YLabel.String = sprintf('%s%s',hObject.YLabel.String,' [RANKS]');
    end
    hObject.UserData.doRank = ~hObject.UserData.doRank;
end
%-------------------------------------------------------------------------------

end

%-------------------------------------------------------------------------------
function Dj = TargetDistances(X,targetInd,distMetric)
    % Distances from row targetInd of X to every row of X, defined as in BF_pdist
    % (including the use of overlapping good values when there are NaNs, and the
    % maximum distance for an undefined correlation), without the full pairwise matrix
    switch distMetric
    case {'euclidean','Euclidean'}
        d = X - X(targetInd,:);
        isGood = ~isnan(d);
        d(~isGood) = 0;
        % Rescale for the number of good values, by sqrt(numFeatures/numGood)
        Dj = sqrt(sum(d.^2,2)*size(X,2)./sum(isGood,2));
    case {'corr_fast','abscorr_fast'}
        if any(isnan(X(:)))
            % (BF_NaNCov approximates the correlations when there are NaNs, and only
            % has a full-matrix form)
            Dfull = BF_pdist(X,distMetric,false,[],true);
            Dj = Dfull(:,targetInd);
        else
            Dj = TargetDistances(X,targetInd,'Pearson');
            if strcmp(distMetric,'abscorr_fast')
                Dj = 1 - abs(1 - Dj);
            end
        end
    case {'Spearman','abscorr','Pearson','PearsonCorr','corr','correlation'}
        if ismember(distMetric,{'Spearman','abscorr'})
            corrType = 'Spearman';
        else
            corrType = 'Pearson';
        end
        if any(isnan(X(:)))
            rowsHow = 'pairwise';
        else
            rowsHow = 'all';
        end
        rho = corr(X(targetInd,:)',X','Type',corrType,'Rows',rowsHow)';
        if strcmp(distMetric,'abscorr')
            Dj = 1 - abs(rho);
            Dj(isnan(rho)) = 1;
        else
            Dj = 1 - rho;
            Dj(isnan(rho)) = 2;
        end
        Dj(targetInd) = 0;
    otherwise
        error('Unknown distance metric: %s',distMetric);
    end
end
