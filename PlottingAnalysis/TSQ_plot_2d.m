% ------------------------------------------------------------------------------
% TSQ_plot_2d
% ------------------------------------------------------------------------------
% 
% Plots the dataset in a two-dimensional space
% e.g., that of two chosen operations, or two principal components.
% 
%---HISTORY
% Borrows from TSQ_pca plotting routines
% Ben Fulcher 24/3/2010
% Ben Fulcher 28/4/2010: added F input
% Ben Fulcher 19/10/2010: added TrainTest input: cell to plot seperately
%                TrainTest = {traini,testi}; if just one component then a
%                subset
% Ben Fulcher 20/10/2010: added annotatep option: number of time series
%                annotations to make (default = 0)
% Ben Fulcher 20/10/2010: also added keepksdensities option
% Ben Fulcher 27/10/2010: added lossmeth option (choose how to calculate
%                           loss)
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

function TSQ_plot_2d(Features,DataInfo,TrainTest,annotatep,keepksdensities,lossmeth)

% Features should be a Nx2 vector of where to plot each of the N data objects in the two-dimensional space
if nargin < 1
    error('You must provide two-dimensional feature vectors for the data.')
end

% DataInfo should be a structure array with all the information about the data (same length as Features)
% Group should be a field in this structure array

if nargin < 3 || isempty(TrainTest)
    TrainTest = {};
end

if nargin < 4 || isempty(annotatep)
    annotatep = struct('n',0);
end

if nargin < 5 || isempty(keepksdensities)
    keepksdensities = 1;
end

if nargin < 6 || isempty(lossmeth)
    lossmeth = 'linclass';
end

MakeFigure = 1; % default is to plot on a brand new figure('color','w')

% ------------------------------------------------------------------------------
%% Load data
% ------------------------------------------------------------------------------
% Data is not loaded, now it must be provided

labels = DataInfo.labels; % Feature labels
if isstruct(annotatep) || annotatep > 0
    DataLabels = DataInfo.DataLabels; % We need data labels
end
GroupNames = DataInfo.GroupNames;
GroupIndices = DataInfo.GroupIndices;
TimeSeriesData = DataInfo.TimeSeriesData;
NumGroups = length(GroupNames);

% ------------------------------------------------------------------------------
%% Subset
% ------------------------------------------------------------------------------
% Only use a subset of the full matrix
if (length(TrainTest)==1 || ~iscell(TrainTest))
    if iscell(TrainTest)
        rss = TrainTest{1}; % row subset
    else
        rss = TrainTest;
    end
    fprintf(1,'Subset rows from %u to %u.\n',size(Features,1),length(rss))
    Features = Features(rss,:);
    % GroupIndices refers to indicies of the full matrix, we want to convert to subset
    % matrix
    % for each subset index, label with it's group number
    grpnum = zeros(length(rss),1);
    for i = 1:length(rss)
        grpnum(i) = find(cellfun(@(x)ismember(rss(i),x),GroupIndices));
    end
    % now we make a new GroupIndices
    ugrpnum = unique(grpnum);
    GroupIndices = cell(length(ugrpnum),1);
    for i = 1:length(GroupIndices)
        GroupIndices{i} = find(grpnum==ugrpnum(i));
    end
    TrainTest = []; % make empty so don't plot TrainTest groups later
end

% ------------------------------------------------------------------------------
% Compute classification rates
% ------------------------------------------------------------------------------
loss = zeros(3,1); % loss1, loss2, lossboth
gig = BF_ToGroup(GroupIndices); % Convert GroupIndices to group form
switch lossmeth
    case 'linclass'
        kfold = 10;
        nrepeats = 5;
        Classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
                        sum(ytest == classify(Xtest,XTrain,yTrain,'linear'))/length(ytest);
        try
            loss(1) = mean(Classify_fn(Features(:,1),gig,Features(:,1),gig));
            loss(2) = mean(Classify_fn(Features(:,2),gig,Features(:,2),gig));
            loss(3) = mean(Classify_fn(Features(:,1:2),gig,Features(:,1:2),gig));
        catch emsg
            fprintf(1,'%s\n',emsg.message);
        end
        fprintf(1,'Linear in-sample classification rates computed\n');
    case {'knn','knn_matlab'}
        k = 3;
        nrepeats = 10;
        try
            lf1 = TSQ_cfnerr('knn',k,Features(:,1),gig,[],{'kfold',10,nrepeats});
            loss(1,:) = [mean(lf1),std(lf1)];
            lf2 = TSQ_cfnerr('knn',k,Features(:,2),gig,[],{'kfold',10,nrepeats});
            loss(2,:) = [mean(lf2),std(lf2)];
            lf3 = TSQ_cfnerr('knn',k,Features(:,1:2),gig,[],{'kfold',10,nrepeats});
            loss(3,:) = [mean(lf3),std(lf3)];
        catch emsg
            fprintf(1,'%s\n',emsg);
        end
end

% ------------------------------------------------------------------------------
%% Plot
% ------------------------------------------------------------------------------
if MakeFigure % can set extras.MakeFigure = 0 to plot within a given setting
    figure('color','w'); box('on'); % white figure
end

% Set colors
if isstruct(annotatep) && isfield(annotatep,'cmap')
    if ischar(annotatep.cmap)
        c = BF_getcmap(annotatep.cmap,NumGroups,1);
    else
        c = annotatep.cmap; % specify the cell itself
    end
else
    if NumGroups < 10
        c = BF_getcmap('set1',NumGroups,1);
    elseif NumGroups <= 12
        c = BF_getcmap('set3',NumGroups,1);
    elseif NumGroups<=22
        c = [BF_getcmap('set1',NumGroups,1); ...
                    BF_getcmap('set3',NumGroups,1)];
    elseif NumGroups <= 50
        c = mat2cell(jet(NumGroups),ones(NumGroups,1));
    else
        error('There aren''t enough colors in the rainbow to plot this many classes!')
    end
end
if (NumGroups == 1)
    c = {'k'}; % Just use black...
end

% ------------------------------------------------------------------------------
%% Plot distributions
% ------------------------------------------------------------------------------
if keepksdensities
    subplot(4,4,1:3); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:NumGroups
        fr = plot_ks(Features(GroupIndices{i},1),c{i},0);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ylim',[minn,maxx]);
    
    subplot(4,4,[8,12,16]); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:NumGroups
        fr = plot_ks(Features(GroupIndices{i},2),c{i},1);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'xlim',[minn,maxx]);
end

% ------------------------------------------------------------------------------
%% Set up a 2D plot
% ------------------------------------------------------------------------------
if keepksdensities
    subplot(4,4,[5:7,9:11,13:15]); box('on');
end
hold on;

if isfield(annotatep,'TheMarkerSize');
    TheMarkerSize = annotatep.TheMarkerSize; % specify custom marker size
else
    if isempty(TrainTest)
        TheMarkerSize = 12; % Marker size for '.'
    else
        TheMarkerSize = 5; % Marker size for 'o' and 's'
    end
end

if isempty(TrainTest)
    for i = 1:NumGroups
        plot(Features(GroupIndices{i},1),Features(GroupIndices{i},2),'.','color',c{i},'MarkerSize',TheMarkerSize)
    end
else % Plot training and test data differently
    for j = 1:length(TrainTest)
        for i = 1:NumGroups
            if (j==1)
                % Training data
                plot(Features(intersect(GroupIndices{i},TrainTest{j}),1),Features(intersect(GroupIndices{i},TrainTest{j}),2),...
                        'ok','MarkerFaceColor',c{i},'MarkerSize',TheMarkerSize)
            else
                % Test data
                plot(Features(intersect(GroupIndices{i},TrainTest{j}),1),Features(intersect(GroupIndices{i},TrainTest{j}),2),...
                        'sk','MarkerFaceColor',c{i},'MarkerSize',TheMarkerSize)
            end
        end
    end
end

%% Plot cluster centres
% for i = 1:NumGroups
%     cc = median(F(gi{i},:));
%     plot(cc(1),cc(2),'o','color',c{i},'MarkerFaceColor',c{i}+(1-c{i})*0.5,...
%                     'MarkerSize',10,'LineWidth',2);
% end

% ------------------------------------------------------------------------------
%% Plot a classify boundary?
% ------------------------------------------------------------------------------
if (NumGroups == 2) && strcmp(lossmeth,'linclass');
    modeorder = 'linear'; % or 'quadratic'
    
    xlim = get(gca,'XLim'); ylim=get(gca,'YLim');
    group = BF_ToGroup(GroupIndices);
    [X, Y] = meshgrid(linspace(xlim(1),xlim(2),200),linspace(ylim(1),ylim(2),200));
    X = X(:); Y = Y(:);
    [~,~,~,~,coeff] = classify([X Y],Features(:,1:2), group, modeorder);
    
    hold on;
    K = coeff(1,2).const; L = coeff(1,2).linear;
    if strcmp(modeorder,'linear')
        Q = zeros(2,2);
    else
        Q = coeff(1,2).quadratic;
    end
    f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
    h2 = ezplot(f,[xlim(1), xlim(2), ylim(1), ylim(2)]);
    set(h2,'LineStyle','--','color','k','LineWidth',2)
    % ezplot overrides the figure title. Reinstate:
%     title(['combined loss (' lossmeth ') = ' num2str(round(loss(3)*100)) '%']);
%     xlabel([labels{mr(1)} ' [' keywords{mr(1)} '] -- loss = ' num2str(round(loss(1)*100)) '%'],'interpreter','none')
%     ylabel('')
end

% ------------------------------------------------------------------------------
%% Label Axes
% ------------------------------------------------------------------------------
title(sprintf('Combined misclassification rate (%s) = %.2f%%',lossmeth, ...
                    round(loss(3,1)*100)),'interpreter','none');

LabelText = cell(2,1);
for i = 1:2
    LabelText{i} = sprintf('%s (%.2f %%)',labels{i}, ...
                            round(loss(i,1)*100)); %,round(loss(i,2)*100));
end
xlabel(LabelText{1},'interpreter','none')
ylabel(LabelText{2},'interpreter','none')

% Set Legend
if isempty(TrainTest)
    legs = cell(NumGroups,1);
    for i = 1:NumGroups
        legs{i} = sprintf('%s (%u)',GroupNames{i},length(GroupIndices{i}));
    end
else
    legs = cell(NumGroups*2,1);
    for i = 1:NumGroups
        legs{i} = sprintf('%s train (%u)',GroupNames{i},length(intersect(GroupIndices{i},TrainTest{1})));
        legs{NumGroups+i} = sprintf('%s test (%u)',GroupNames{i},length(intersect(GroupIndices{i},TrainTest{2})));
    end
end
legend(legs);

% ------------------------------------------------------------------------------
%% Annotate time-series data
% ------------------------------------------------------------------------------
if isempty(TimeSeriesData)
    % Only attempt to annotate if time-series data is provided
    return
end
% Set parameters
if isfield(annotatep,'maxL')
    maxL = annotatep.maxL;
else
    maxL = 300; % length of annotated time series segments
end
if isfield(annotatep,'UserInput')
    UserInput = annotatep.UserInput;
else
    UserInput = 1; % user input points rather than randomly chosen
end
if isfield(annotatep,'fdim')
    fdim = annotatep.fdim;
else
    fdim = [0.30,0.08]; % width, height
end
if isfield(annotatep,'TextAnnotation')
    TextAnnotation = annotatep.TextAnnotation; % 'filename','tsid','none'
else
    TextAnnotation = 'none'; % no annotations by default
end
if isfield(annotatep,'TheLineWidth')
    TheLineWidth = annotatep.TheLineWidth;
else
    TheLineWidth = 0.8; % % line width for time series
end
NumAnnotations = annotatep.n;

pxlim = get(gca,'xlim'); % plot limits
pylim = get(gca,'ylim'); % plot limits
pwidth = diff(pxlim); % plot width
pheight = diff(pylim); % plot height
AlreadyPicked = zeros(NumAnnotations,2); % record those already picked
PlotCircle = 1; % magenta circle around annotated points

% produce xy points
xy = cell(NumGroups,1);
for i = 1:NumGroups
    xy{i} = [Features(GroupIndices{i},1),Features(GroupIndices{i},2)];
end

if ~UserInput % points to annotate are randomly picked
    if NumAnnotations == length(DataLabels) % annotate all
        fprintf(1,'Annotate all\n')
        for j = 1:NumAnnotations
            TheGroup = find(cellfun(@(x)ismember(j,x),GroupIndices));
            AlreadyPicked(j,1) = TheGroup;
            AlreadyPicked(j,2) = find(GroupIndices{TheGroup}==j);
        end
    else
        AlreadyPicked(:,1) = round(linspace(1,NumGroups,NumAnnotations));
        randperms = cellfun(@(x)randperm(length(x)),GroupIndices,'UniformOutput',0);
        counters = ones(NumGroups,1);
        for j = 1:NumAnnotations
            AlreadyPicked(j,2) = randperms{AlreadyPicked(j,1)}(counters(AlreadyPicked(j,1))); % random element of the group
            counters(AlreadyPicked(j,1)) = counters(AlreadyPicked(j,1))+1;
        end
    end
end

for j = 1:NumAnnotations
    if UserInput
        point = ginput(1);
        iplot = ClosestPoint_ginput(xy,point); % find closest actual point to input point
        TheGroup = iplot(1); % want this group
        itsme = iplot(2); % and this index
        AlreadyPicked(j,:) = [TheGroup,itsme];
    else
        TheGroup = AlreadyPicked(j,1);
        itsme = AlreadyPicked(j,2);
    end
    
    if (j > 1) && any(sum(abs(AlreadyPicked(1:j-1,:)-repmat(AlreadyPicked(j,:),j-1,1)),2)==0)
        % Same one has already been picked, don't plot it again
        continue
    end
    
    PlotPoint = xy{TheGroup}(itsme,:);
    % fn = DataLabels{GroupIndices{TheGroup}(itsme)}; % filename of timeseries to plot
    % ts = dlmread(fn);
    ts = TimeSeriesData{GroupIndices{TheGroup}(itsme)}; % filename of timeseries to plot
    if ~isempty(maxL)
        ts = ts(1:min(maxL,end));
    end
    
    if PlotCircle
        plot(PlotPoint(1),PlotPoint(2),'om'); % Plot a magenta circle around the target point
    end
    
    switch TextAnnotation
    case 'filename'
        % Annotate text with names of datapoints:
        text(PlotPoint(1),PlotPoint(2)-0.01*pheight,fn,'interpreter','none','FontSize',8);
    case 'tsid'
        % Annotate text with ts_id:
        text(PlotPoint(1),PlotPoint(2)-0.01*pheight,num2str(ts_ids_keep(GroupIndices{TheGroup}(itsme))),'interpreter','none','FontSize',8);
    end
    
    % Adjust if annotation goes off axis x-limits
    px = PlotPoint(1)+[-fdim(1)*pwidth/2,+fdim(1)*pwidth/2];
    if px(1) < pxlim(1), px(1) = pxlim(1); end % can't plot off left side of plot
    if px(2) > pxlim(2), px(1) = pxlim(2)-fdim(1)*pwidth; end % can't plot off right side of plot
    
    % Adjust if annotation goes above maximum y-limits
    py = PlotPoint(2)+[0,fdim(2)*pheight];
    if py(2) > pylim(2)
        py(1) = pylim(2)-fdim(2)*pheight;
    end
    
    % Annotate the time series
    plot(px(1)+linspace(0,fdim(1)*pwidth,length(ts)),...
            py(1)+fdim(2)*pheight*(ts-min(ts))/(max(ts)-min(ts)),...
                '-','color',c{TheGroup},'LineWidth',TheLineWidth);
end

function [fr, xr] = plot_ks(v,c,swap)
    % Vector v is the vector of a given group
    % c is the color
    [f, x] = ksdensity(v,linspace(min(v),max(v),1000),'function','pdf');
    r = zeros(length(v),1);
    for m = 1:length(v); r(m)=find(x>=v(m),1,'first'); end
    r = sort(r,'ascend');
    fr = f(r); xr = x(r);
    if swap
        plot(f,x,'color',c);
        plot(fr,xr,'.','color',c)
    else
        plot(x,f,'color',c);
        plot(xr,fr,'.','color',c)
    end
end
	
end