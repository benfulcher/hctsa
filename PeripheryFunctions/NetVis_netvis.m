function NetVis_netvis(A,varargin)
% NetVis_netvis     Network visualization.
%
% Adjacency matrix to be plotted is A. A range of other options are available
% to set visualization and annotation settings.
%
% If not a single connected component, set k > 0, see below.
%
% Assumptions:
% 1. Adjacency matrix is symmetric.
% 2. Adjacency matrix is not sparse.
%
%---EXAMPLE USAGE:
% NetVis_netvis(A,'k',0.01,...
%                'textLabels',textLabels,...
%                'linkThresh',[0.9,0.8,0.7,0.6],...
%                'nodeLabels',nodeLabels,...
%                'dataLabels',dataLabels,...
%                'colorMap','set1');

% ------------------------------------------------------------------------------
% Function designed by Ben Fulcher around the mechanics of the visualization
% algorithm by Dr. Dave Smith, 28/5/2011
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------

% Adjacency matrix, A:
if nargin < 1
    error('Input a valid adjacency matrix');
end
numNodes = length(A);

% Check remaining inputs using an inputParser:
inputP = inputParser;

% Parameter k
% Typically k = 0.01 is a largish value - for fully connected
% graphs k = 0 is valid, k = 1/numNodes is more usual
default_k = 0;
check_k = @(x) validateattributes(x,{'numeric'},{'positive'});
addOptional(inputP,'k',default_k,check_k);

% Text labels for each node, textLabels
default_textLabels = {};
addOptional(inputP,'textLabels',default_textLabels,@iscell);

% linkThresholds, linkThresh (for plotting how-strong links)
default_linkThresh = [];
addOptional(inputP,'linkThresh',default_linkThresh,@isnumeric);

% nodeLabels, A group label (number) for each node
default_nodeLabels = ones(numNodes,1);
addOptional(inputP,'nodeLabels',default_nodeLabels,@isnumeric);

% labelLength
default_labelLength = 0;
check_labelLength = @(x) validateattributes(x,{'numeric'},{'positive'});
addOptional(inputP,'labelLength',default_labelLength,check_labelLength);

% dataLabels, A piece of time-series data for each node
default_dataLabels = cell(numNodes,1);
addOptional(inputP,'dataLabels',default_dataLabels,@iscell);

% Number of iterations, repeats: nits
default_nits = [2000,1]; % [maxIter,numRepeats];
check_nits = @(x) isnumeric(x) && length(x)==2;
addOptional(inputP,'nits',default_nits,check_nits);

% ---PARAMETERS:
% Color map, colorMap
default_cmap = 'set2';
addParameter(inputP,'colorMap',default_cmap,@ischar);

% Tag style, tagStyle, how to display text labels
default_tagStyle = 'k'; % 'k' or 'boxedtext'
addParameter(inputP,'tagStyle',default_tagStyle,@ischar);

% Threshold by proportion or actual threshold on the weights?
default_pthresh = []; % set default below, after par_OperationCodeString
addParameter(inputP,'thresholdByProportion',default_pthresh,@isnumeric);

% Make a new figure?:
default_makeFigure = true;
addParameter(inputP,'makeFigure',default_makeFigure,@islogical);

% Marker size, can specify as a scalar, or as a vector with a custom value for
% each node
default_nodeSize = 12;
addParameter(inputP,'nodeSize',default_nodeSize,@isnumeric);

%-------------------------------------------------------------------------------
%% Parse inputs:
parse(inputP,varargin{:});

% Convert back to variables:
k = inputP.Results.k;
textLabels = inputP.Results.textLabels;
linkThresh = inputP.Results.linkThresh;
nodeLabels = inputP.Results.nodeLabels;
labelLength = inputP.Results.labelLength;
dataLabels = inputP.Results.dataLabels;
nits = inputP.Results.nits;
colorMap = inputP.Results.colorMap;
tagStyle = inputP.Results.tagStyle;
thresholdByProportion = inputP.Results.thresholdByProportion;
makeFigure = inputP.Results.makeFigure;
nodeSize = inputP.Results.nodeSize;

numGroups = length(unique(nodeLabels)); % Number of different (colored) groups to plot

% Set default of thresholdByProportion depending on the link thresholds:
if isempty(thresholdByProportion)
    if any(linkThresh > 1)
        thresholdByProportion = 0;
        fprintf(1,'Thresholding by absolute\n');
    else
        thresholdByProportion = 1;
        fprintf(1,'Thresholding by proportion\n');
    end
end

fprintf(1,'Visualizing a network with %u nodes\n',numNodes);

% ------------------------------------------------------------------------------

fprintf(1,'--------Generating a network visualization with k = %u--------\n',k);

% ------------------------------------------------------------------------------
%% Determine a distance similarity matrix A
% ------------------------------------------------------------------------------
if isempty(linkThresh)
    % Don't threshold -- plot all links
    Ath = cell(1);
    Ath{1} = (A > 0);
    linkThresh = 1;
else
    Ath = cell(length(linkThresh),1);
    for i = 1:length(linkThresh);
        if thresholdByProportion
            if i==1
                Ath{i} = (A >= quantile(A(:),linkThresh(i))); % threshold
            else
                Ath{i} = (A >= quantile(A(:),linkThresh(i)) & A < quantile(A(:),linkThresh(i-1))); % threshold
            end
        else % threshold by absolute
            if i == 1
                Ath{i} = (A >= linkThresh(i)); % threshold
            else
                Ath{i} = (A >= linkThresh(i) & A < linkThresh(i-1)); % threshold
            end
        end
    %     Ath = (A>th1 & A<=th2); % threshold (2) % absolute threshold
    end
end

% Find versions that contain links:
anyLinks = find(cellfun(@(x)any(x(:)),Ath));
% Ath = Ath(anyLinks);


% ------------------------------------------------------------------------------
% Prepare colors
% ------------------------------------------------------------------------------

if iscell(colorMap)
    c = colorMap; % specified a cell of colors
else % specified the name of a color set to load
    c = BF_getcmap(colorMap,numGroups,1);
end
if length(c) < numGroups
    fprintf(1,'There aren''t enough colors to plot %u groups of the data\n',numGroups);
    fprintf(1,'Changing to jet colormap\n');
    myJet = jet(numGroups);
    c = cell(numGroups,1);
    for i = 1:numGroups
        c{i} = myJet(i,:);
    end
end

% Color of links
maxBrighten = 0.8;
clinks = BF_getcmap('blues',3,1);
clinks = clinks{3};
clinks = arrayfun(@(x) brighten(clinks,x),linspace(0,maxBrighten,length(linkThresh)),'UniformOutput',0);
% clinks = flipud(clinks); % dark colours first
% clinks = BF_getcmap('blues',max(3,length(linkThresh)),1);
% clinks = flipud(clinks); % dark colours first

% Text color:
tc = cell(length(c),1);
for i = 1:length(c)
    if sum(c{i}) < 2 % Dark background color:
        tc{i} = [1,1,1]; % White text
    else
        tc{i} = [0,0,0]; % Black text
    end
end

% ------------------------------------------------------------------------------
% Make short text labels
% ------------------------------------------------------------------------------
nodeText = cell(numNodes,1);
if ~isempty(textLabels)
    if labelLength > 0 % crop labels to a maximum length
        for i = 1:numNodes
           nodeText{i} = textLabels{i}(1:min(labelLength,end));
        end
    else % keep full label lengths
        nodeText = textLabels;
    end
end

% ------------------------------------------------------------------------------
% Set optimization parameters:
% ------------------------------------------------------------------------------
maxIter = nits(1); % iterations per run;
numRepeats = nits(2); % number of repeats
numApprox = 15; % Number of approximations used in the evaluation of the repulsion forces

% ------------------------------------------------------------------------------
% Deduct the diagonal entries
% ------------------------------------------------------------------------------
B = A - diag(diag(A));
L = diag(sum(B,2)) - B; % Laplacian
E = L + k*numNodes*speye(numNodes);
R = chol(sparse(E));
Rtr = R';

Energy = zeros(numRepeats,3); % energies, gradients
xy = cell(numRepeats,1); % sets of node positions (x,y)

fprintf(1,'Computing network visualization over %u iterations and %u repeats\n',maxIter,numRepeats);

for jman = 1:numRepeats
    % Random intial values - or could seed with results from an earlier layout
    x = rand(numNodes,1); x = x - mean(x);
    y = rand(numNodes,1); y = y - mean(y);

    % Run through the mechanics of the visualization method:
    for iter = 1:maxIter  % Could use a threshold on the gradient here
        [f,g] = NetVis_Jtilda2d(x,y,numApprox); % Gets the repulsion forces using an approximation
        J = [(E*x-f);(E*y-g)]; % This is the full gradient (attraction - repulsion)
        % This is only to plot the residual (so can be removed for speed)
        x = R\(Rtr\f); x = x-mean(x); % Update and center
        y = R\(Rtr\g); y = y-mean(y); % Update and center
    end

    % Give information about convergence:
    DD = BF_pdist([x,y],'euclidean');
    Energy(jman,1) = sum(sum((A+k).*DD.^2 - 2*DD));
    Energy(jman,2) = log(J'*J);
    RR = corrcoef(DD,A);
    Energy(jman,3) = RR(2,1);
    xy{jman} = [x,y];

    if numRepeats > 1
        repeatText = sprintf('Repeat %u/%u: ',jman,numRepeats);
    else
        repeatText = '';
    end

    fprintf(1,['%sE = %.2g [%.2g] (grad = %.2g)\n' ...
                    'corr(input_distance,2d_euclidean_distance) = %.2g\n'],...
                repeatText, ...
                Energy(jman,1), ...
                Energy(jman,1) - min(Energy(1:jman-1,1)), ...
                Energy(jman,2), ...
                Energy(jman,3));
end

% Now run the best one further (could do this starting with seeds from
% previous layout) but I won't do it.
[~,theBest] = min(Energy(:,1));
x = xy{theBest}(:,1);
y = xy{theBest}(:,2);

% ------------------------------------------------------------------------------
%% PLOT IT!!
% ------------------------------------------------------------------------------
% Make a new figure?
if makeFigure
    figure('Position',[10,150,500,500],'color','w');
end

% Aesthetics
ax = gca;
ax.XTick = [];
ax.XTickLabel = {};
ax.YTick = [];
ax.YTickLabel = {};

% ------------------------------------------------------------------------------
% PLOT LINKS
% ------------------------------------------------------------------------------
hold off
% Go through each thresholded adjacency matrix (mutually exclusive links), Ath
% and plot the links in them a different color
sizeDecline = linspace(1,0.5,length(linkThresh));
ordering = sort(anyLinks,2,'descend');
handles = cell(length(ordering),1);
for i = 1:length(ordering) % length(Ath):-1:1 % plot the weakest links first
    ii = ordering(i);
    [X,Y] = gplot(Ath{ii},[x y]); %,'color',clinks{2}); %,'color',c{2}); % links
    handles{i} = plot(X(:),Y(:),'Color',clinks{ii},'LineWidth',2*sizeDecline(ii)); % (length(A)-i+1)*0.8
    hold on
end

% ------------------------------------------------------------------------------
% PLOT NODES
% ------------------------------------------------------------------------------
pwidth = diff(get(gca,'XLim')); % plot width
pheight = diff(get(gca,'YLim')); % plot height
fdim = [0.2,0.05]; % fractional width/height of time-series annotations
maxL = 100; % maximum number of annotated samples
for i = 1:numGroups
    if length(nodeSize) > 1
        % sizes specified for each node as a vector
        Nodes_in_i = find(nodeLabels==i);
        for j = 1:length(Nodes_in_i)
            % plot nodes individually
            plot(x(Nodes_in_i(j)),y(Nodes_in_i(j)),'o','MarkerFaceColor',c{i},...
                'MarkerEdgeColor','k','markerSize',nodeSize(Nodes_in_i(j)));
        end

    else % default Marker size = 7
        % Plot all together
        plot(x(nodeLabels==i),y(nodeLabels==i),'o', ...
                'MarkerFaceColor',c{i},'MarkerEdgeColor','k', ...
                'LineWidth',1.5,'markerSize',nodeSize); %,'MarkerEdgeColor','k'); % nodes
    end

    % Annotate time series to the plot:
    if ~isempty(dataLabels)
        % do it as a proportion (fdim) of plot size
        Nodes_in_i = find(nodeLabels==i);
        for j = 1:length(Nodes_in_i)
            ts = dataLabels{Nodes_in_i(j)}(1:min(maxL,end));
            plot(x(Nodes_in_i(j)) + linspace(0,pwidth*fdim(1),length(ts)),...
                    y(Nodes_in_i(j)) - pheight*fdim(2) + pheight*fdim(2)*(ts-min(ts))/(max(ts)-min(ts)),...
                    'Color',c{i},'LineWidth',1.0)
        end
    end

    % ------------------------------------------------------------------------------
    % Text annotations:
    % ------------------------------------------------------------------------------
    inc = 0.0; textFontSize = 11;
    switch tagStyle
    case 'coloredtext' % colored text
        for j = find(nodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                 'interpreter','none','Color',c{i});
        end
    case 'boxedtext' % boxed, colored text
        for j = find(nodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                    'interpreter','none','BackGroundColor',c{i},'Color',tc{i});
        end
    case 'k' % black text labels
        for j = find(nodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                 'interpreter','none','Color','k');
        end
    end
end

%-------------------------------------------------------------------------------
% Set a legend for links:
legend([handles{:}],arrayfun(@(x)num2str(linkThresh(x)),sort(anyLinks,'descend'),'UniformOutput',false))

% ------------------------------------------------------------------------------
% Tidy up:
% ------------------------------------------------------------------------------
axis square
set(gca,'XTick',[],'XTickLabel',{},'YTick',[],'YTickLabel',{});
box('on')

end
