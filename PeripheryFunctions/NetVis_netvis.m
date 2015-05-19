% Network Visualization
% 
% Adjacency matrix to be plotted is A.
% 
% If not a single connected component, set k > 0, see below.
% 
% Assumptions:
% 1. Adjacency matrix is symmetric.
% 2. Adjacency matrix is not sparse.
% 
% ------------------------------------------------------------------------------
% Mechanics of network code from Dr. Dave Smith, 28/5/2011
% Alterations made by Ben Fulcher
% ------------------------------------------------------------------------------

function BNetVis(A,k,labelscl,linkThresh,NodeLabels,labelLength,nits,opts)


% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 1
    error('Need an adjacency matrix for input...');
end
numNodes = length(A);
fprintf(1,'Visualizing a network with %u nodes\n');

% Parameter k %typically k = 0.01 is a largish value - for fully connected
% graphs k = 0 is valid, k = 1/numNodes is more usual
if nargin < 2 || isempty(k)
    k = 0;
end

if nargin < 3
    labelscl = {}; % Text labels for each node
end

if nargin < 4 || isempty(linkThresh)
    linkThresh = []; % [0.1,0.2]; % 10% of links are links, others dotted
end

if nargin < 5 || isempty(NodeLabels)
    NodeLabels = ones(numNodes,1); % {(1:numNodes)};
    % A label (number) for each node
end
NumGroups = length(unique(NodeLabels)); % Number of different (colored) groups to plot

if nargin < 6 || isempty(labelLength)
    labelLength = 0;
end

if nargin < 7 || isempty(nits)
    nits = [1000,1]; % [maxIter,numRepeats]
end

if nargin < 8 || isempty(opts)
    opts = 'set2';
end

if ischar(opts),
    cmap = opts;
elseif isstruct(opts) && isfield(opts,'cmap')
    cmap = opts.cmap;
end

if isstruct(opts) && isfield(opts,'tagStyle')
    tagStyle = opts.tagStyle;
else
    tagStyle = 'k'; % also 'boxedtext' -- how to display text labels
end 

% Threshold by proportion or actual threshold on the weights?
if isstruct(opts) && isfield(opts,'pthresh')
    thresholdByProportion = opts.pthresh;
else
    if any(linkThresh > 1)
        thresholdByProportion = 0;
        fprintf(1,'Thresholding by absolute\n')
    else
        thresholdByProportion = 1;
        fprintf(1,'Thresholding by proportion\n')
    end
end


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

if iscell(cmap)
    c = cmap; % specified a cell of colors
else % specified the name of a color set to load
    c = BF_getcmap(cmap,NumGroups,1);
end
if length(c) < NumGroups
    fprintf(1,'There aren''t enough colors to plot %u groups of the data\n',NumGroups);
    fprintf(1,'Changing to jet colormap\n');
    myJet = jet(NumGroups);
    c = cell(NumGroups,1);
    for i = 1:NumGroups
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
if ~isempty(labelscl)
    if labelLength > 0 % crop labels to a maximum length
        for i = 1:numNodes
           nodeText{i} = labelscl{i}(1:min(labelLength,end));
        end
    else % keep full label lengths
        nodeText = labelscl;
    end
end

% ------------------------------------------------------------------------------
% Number of approximations used in the evaluation of the repulsion forces:
% ------------------------------------------------------------------------------
numApprox = 15;

% ------------------------------------------------------------------------------
% Maximum number of iterations:
% ------------------------------------------------------------------------------
maxIter = nits(1); % 100;
numRepeats = nits(2); % number of repeats

% ------------------------------------------------------------------------------
% Deduct the diagonal entries
% ------------------------------------------------------------------------------
B = A - diag(diag(A));
L = diag(sum(B,2)) - B; % Laplacian
E = L + k*numNodes*speye(numNodes);
[R, p] = chol(sparse(E));
Rtr  = R';

Energy = zeros(numRepeats,3); % energies, gradients
xy = cell(numRepeats,1); % sets of node positions (x,y)

fprintf(1,'Computing network visualization over %u iterations and %u repeats\n',maxIter,numRepeats);

for jman = 1:numRepeats
    % Random intial values - or could seed with results from an earlier layout
    x = rand(numNodes,1); x = x - mean(x);
    y = rand(numNodes,1); y = y - mean(y);
    
    % Run through the mechanics of the visualization method:
    for iter = 1:maxIter  % Could use a threshold on the gradient here
        [f g] = NetVis_Jtilda2d(x,y,numApprox); % Gets the repulsion forces using an approximation
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
    
    fprintf(1,['Iteration %u/%u: E = %.2g [%.2g] (grad = %.2g),' ...
                    ' corr(input_distance,2d_euclidean_distance) = %.2g)\n'],...
                jman,numRepeats, ...
                Energy(jman,1), ...
                Energy(jman,1) - min(Energy(1:jman-1,1)), ...
                Energy(jman,2), ...
                Energy(jman,3));
end

% Now run the best one further (could do this starting with seeds from
% previous layout) but I won't do it.
[minE,theBest] = min(Energy(:,1));
x = xy{theBest}(:,1);
y = xy{theBest}(:,2);

% ------------------------------------------------------------------------------
%% PLOT IT!!
% ------------------------------------------------------------------------------
% Make a new figure?
if isfield(opts,'makeFigure')
    makeFigure = opts.makeFigure;
else
    makeFigure = 1;
end

if makeFigure
    figure('Position',[10,150,500,500],'color','w');
end

% Aesthetics
set(gca,'XTick',[],'XTickLabel',{},'YTick',[],'YTickLabel',{});
% if ~isempty(tsx)
%     fdim = [0.25,0.06]; % width, height
% end

% ------------------------------------------------------------------------------
% PLOT LINKS
% ------------------------------------------------------------------------------
hold off
% Go through each thresholded adjacency matrix (mutually exclusive links), Ath
% and plot the links in them a different color
sizeDecline = linspace(1,0.5,length(linkThresh));
for i = sort(anyLinks','descend') % length(Ath):-1:1 % plot the weakest links first
    try
        [X,Y] = gplot(Ath{i},[x y]); %,'color',clinks{2}); %,'color',c{2}); % links
        plot(X(:),Y(:),'Color',clinks{i},'LineWidth',2*sizeDecline(i)); % (length(A)-i+1)*0.8
        hold on
    catch
        keyboard
    end
end
% Set a legend:
legend(arrayfun(@(x)num2str(linkThresh(x)),sort(anyLinks,'descend'),'UniformOutput',0))

% hold off;
% [X,Y] = gplot(Ath{1},[x,y]); %,'color',clinks{2}); %,'color',c{2}); % links
% plot(X(:),Y(:),'Color',clinks{1},'LineWidth',2*0.5^i); % (length(A)-i+1)*0.8
% for i = 1:NumGroups
%     for j = i:NumGroups
%         if i==j
%             TheColor = c{i};
%         else
%             TheColor = clinks{i};
%         end
%         [X,Y] = gplot(Ath{1},[x,y]); %,'color',clinks{2}); %,'color',c{2}); % links
%         plot(X(:),Y(:),'Color',clinks{i},'LineWidth',2*0.5^i); % (length(A)-i+1)*0.8
%     end
% end
% hold on


%     [X,Y] = gplot(A,[x y],'-'); %,'color',c{2}); % links
%     plot(X(:),Y(:),'Color',clinks{1});
%     hold on

% ------------------------------------------------------------------------------
% PLOT NODES
% ------------------------------------------------------------------------------
pwidth = diff(get(gca,'XLim')); % plot width
pheight = diff(get(gca,'YLim')); % plot height
for i = 1:NumGroups
    if isfield(opts,'size') && length(opts.size) > 1
        % sizes specified for each node as a vector
        Nodes_in_i = find(NodeLabels==i);
        for j = 1:length(Nodes_in_i)
            % plot nodes individually
            plot(x(Nodes_in_i(j)),y(Nodes_in_i(j)),'o','MarkerFaceColor',c{i},...
                'MarkerEdgeColor','k','MarkerSize',opts.size(Nodes_in_i(j)));
        end
        
    else % default Marker size = 7
        if isfield(opts,'size')
            msize = opts.size;
        else
            msize = 12;
        end

        % Plot all together
        plot(x(NodeLabels==i),y(NodeLabels==i),'o', ...
                'MarkerFaceColor',c{i},'MarkerEdgeColor','k', ...
                'LineWidth',1.5,'MarkerSize',msize); %,'MarkerEdgeColor','k'); % nodes
    end
    % Annotate time series to the plot:
    % if ~isempty(tsx)
    %     % do it as a proportion (fdim) of plot size
    %     Nodes_in_i = find(NodeLabels==i);
    %     for j = 1:length(Nodes_in_i)
    %         ts = tsx{Nodes_in_i(j)};
    %         plot(x(Nodes_in_i(j)) + linspace(0,pwidth*fdim(1),length(ts)),...
    %                 y(Nodes_in_i(j)) - pheight*fdim(2) + pheight*fdim(2)*(ts-min(ts))/(max(ts)-min(ts)),...
    %                 'Color',c{i},'LineWidth',1.0)
    %     end
    % end
    
    % ------------------------------------------------------------------------------
    % Text annotations:
    % ------------------------------------------------------------------------------
    inc = 0.0; textFontSize = 11;
    switch tagStyle
    case 'coloredtext' % colored text
        for j = find(NodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                 'interpreter','none','Color',c{i});
        end
    case 'boxedtext' % boxed, colored text
        for j = find(NodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                    'interpreter','none','BackGroundColor',c{i},'Color',tc{i});
        end
    case 'k' % black text labels
        for j = find(NodeLabels==i)
            text(x(j)+inc,y(j)+inc,nodeText(j),'FontSize',textFontSize, ...
                 'interpreter','none','Color','k');
        end
    end
end

% ------------------------------------------------------------------------------
% Tidy up:
% ------------------------------------------------------------------------------
axis square
set(gca,'XTick',[],'XTickLabel',{},'YTick',[],'YTickLabel',{});
box('on')

end