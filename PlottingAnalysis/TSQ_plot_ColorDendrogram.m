% TSQ_plot_ColorDendrogram
% 
% Creates a dendrogram for distance matrix d with colours at each leaf
% 
%---INPUTS:
% d is the (NxN) pairwise distance matrix (e.g., computed using pdist, then squareform)
% The Nx1 vector groups gives the index each node belongs to
% The Gx1 cell GroupLabels gives the name of each group (for each of G groups)
% LinkMethod specifies the linkage method ('complete' by default)
% horiz can plot everything horizontally. Don't think this works yet...
% 
%---OUTPUT:
% Generates a formatted and color-labeled dendrogram of your data
% 
%---HISTORY:
% Ben Fulcher 2010, based on code by Dann Fenn
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

function TSQ_plot_ColorDendrogram(d,groups,GroupLabels,LinkMethod,horiz)

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 4 || isempty(LinkMethod) % linkage method
    LinkMethod = 'complete'; % use average clustering by default
end
if nargin < 5 || isempty(horiz)
    horiz = 0; % display horizontally? not really supported yet
end
if horiz==1
    warning('horizontal displays not really supported yet, sorry... :(')
end

% Set up orientation: Orientation
if horiz
    Orientation = 'left';
else
    Orientation = 'top';
end

% Can specify distance matrix as a string: norcl, and will retrieve and
% calculate distances as appropriate...
if ischar(d)
    dmth = 'euclidean';
    
    if strcmp(d,'cl')
        TheFile = 'HCTSA_cl.mat';
    elseif strcmp(d,'norm')
        TheFile = 'HCTSA_N.mat';
    end
    
    load(TheFile,'TS_DataMat');
    F = TS_DataMat; clear TS_DataMat
    
    % Calculate pairwise distances
    d = BF_pdist(F,dmth);
end

% Rescale distances to a maximum of 5, helps visualization
newmaxd = 5;
d = newmaxd * d / max(d);

% Each cell element are indices of that group
if iscell(groups)
    groups = BF_ToGroup(groups);
end

% ------------------------------------------------------------------------------
% Do the linkage:
% ------------------------------------------------------------------------------
Z = linkage(d,LinkMethod);

% ------------------------------------------------------------------------------
% Assign labels
% ------------------------------------------------------------------------------
NumNodes = length(groups); % The number of objects/nodes
NodeLabels = cell(NumNodes,1);
for i=1:NumNodes; NodeLabels{i} = ''; end; % Don't actually label the nodes...

% ------------------------------------------------------------------------------
%% Plot the dendrogram
% ------------------------------------------------------------------------------
figure('color','w'); hold on;
if length(groups)>1000
    [H,~,perm] = dendrogram(Z,0,'Orientation',Orientation,'Labels',NodeLabels);
else % Try optimal leaf order
    fprintf('Running optimal leaf order for %u nodes...',length(groups))
    try
        order = optimalleaforder(Z,d);
        [H,~,perm] = dendrogram(Z,0,'Reorder',order,'Orientation',Orientation,'Labels',NodeLabels);
        fprintf(1,' Done.\n');
    catch
        [H,~,perm] = dendrogram(Z,0,'Orientation',Orientation,'Labels',NodeLabels);
        fprintf(1,' Failed.\n');
    end
end

% Change the line width of the dendrogram
set(H,'LineWidth',0.0002,'Color','k');
set(gca,'FontSize',12);
NumGroups = length(GroupLabels); % number of groups


% ------------------------------------------------------------------------------
% Set up colors
% ------------------------------------------------------------------------------
ng = 6; % use 6 partitions for each color
cmap = colormap(BF_getcmap('blues',ng,0,1));
if NumGroups >= 2
    cmap = [cmap;BF_getcmap('greens',ng,0,1)];
end
if NumGroups >= 3
    cmap = [cmap;BF_getcmap('oranges',ng,0,1)];
end
if NumGroups >= 4
    cmap = [cmap;BF_getcmap('purples',ng,0,1)];
end
if NumGroups >= 5
    cmap = [cmap;BF_getcmap('reds',ng,0,1)];
end
if NumGroups >= 6
    cmap = [cmap;pink(ng)];
end
if NumGroups >= 7
    cmap = [cmap;gray(ng)];
end
if NumGroups >= 8
    cmap = [cmap;BF_getcmap('yelloworangered',ng,0,1)];
end
if NumGroups >= 9
    cmap = [cmap;BF_getcmap('purplebluegreen',ng,0,1)];
end
if NumGroups >= 10
    cmap = [cmap;BF_getcmap('yellowgreenblue',ng,0,1)];
end
if NumGroups >= 11
    error('Too many groups (>10), I can''t handle this!!')
end
colors = cmap;

% ------------------------------------------------------------------------------
% Plot little colored bars to label groups on the dendrogram
% ------------------------------------------------------------------------------
barheight = 0.4;
ugroups = unique(groups);
for u = 1:length(ugroups)
    idx = find(groups==ugroups(u));
    for i = 1:length(idx)
%         h = line([find(perm==idx(i)),find(perm==idx(i))],[-0.02,-0.01],...
%                 'LineWidth',1,'Color',colors(u,:));
%         if i==1, l = [l;h]; end
%         rectangle('Position',[find(perm==idx(i))-0.5,-0.01,1,barheight],.
%         ..
%                 'FaceColor',colors(u,:),'EdgeColor',colors(u,:),'LineWidt
%                 h',0.0001);
        for j = 1:ng
            try
                if horiz
                    rectangle('Position',[-barheight+(j-1)*barheight/ng,find(perm==idx(i))-0.5,barheight/ng,1],...
                        'FaceColor',colors(ng*(u-1)+j,:),'EdgeColor',colors(ng*(u-1)+j,:),'LineWidth',0.0001);
                else
                    rectangle('Position',[find(perm==idx(i))-0.5,-barheight+(j-1)*barheight/ng,1,barheight/ng],...
                        'FaceColor',colors(ng*(u-1)+j,:),'EdgeColor',colors(ng*(u-1)+j,:),'LineWidth',0.0001);
                end
            catch emsg
                keyboard
            end
        end
%         rectangle('Position',[find(perm==idx(i))-0.5,0,1,barheight/6],...
%                 'FaceColor',colors(6*(u-1)+1,:),'EdgeColor',colors(6*(u-1)+1,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+2,:),'EdgeColor',colors(6*(u-1)+2,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,2*barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+3,:),'EdgeColor',colors(6*(u-1)+3,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,3*barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+4,:),'EdgeColor',colors(6*(u-1)+4,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,4*barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+5,:),'EdgeColor',colors(6*(u-1)+5,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,4*barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+5,:),'EdgeColor',colors(6*(u-1)+5,:),'LineWidth',0.0001);
%         rectangle('Position',[find(perm==idx(i))-0.5,5*barheight/6,1,barheight/6],...
%             'FaceColor',colors(6*(u-1)+6,:),'EdgeColor',colors(6*(u-1)+6,:),'LineWidth',0.0001);

    end
end

ylim([-barheight-0.04 max(Z(:,3)) + 0.001]);
ylabel('d','FontSize',20);
hold off;

% ------------------------------------------------------------------------------
%% Write legend, do labeling...
% ------------------------------------------------------------------------------
% legend(GroupLabels)
colors = colors(2:ng:end,:);
y = 0.9;    % vertical height of first colour bar
x = 0.76;   % left corner of each colour bar
w = 0.05;   % width of colour bar
d = 0.045;  % vertical separation between colour bars
for i = 1:length(GroupLabels)
    annotation(gcf,'line',[x x+w],[y-(i-1)*d y-(i-1)*d],'LineWidth',10,'Color',colors(i,:));
    annotation(gcf,'textbox',[x+w+0.005 y-(i-1)*d-0.01 0.1782 0.03361],'String',GroupLabels(i),...
        'FontSize',16,'FitBoxToText','off','LineStyle','none');
end;
set(gcf,'color','w');
