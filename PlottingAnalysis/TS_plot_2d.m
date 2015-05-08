% ------------------------------------------------------------------------------
% TS_plot_2d
% ------------------------------------------------------------------------------
% 
% Plots the dataset in a two-dimensional space
% e.g., that of two chosen operations, or two principal components.
% 
%---INPUTS:
% Features, an Nx2 vector of where to plot each of the N data objects in the
%           two-dimensional space
%
% DataInfo, a structure containing all the information about the data. Fields
%           can include:
%               - labels (feature labels, cols of Features)
%               - DataLabels (labels for each data point, rows of Features)
%               - GroupIndices (a group index for each data point, rows of Features)
%               - GroupNames (name for each group)
%               - TimeSeriesData (cell of vectors containing time-series data)
% 
% trainTest, whether to plot separately indices of training and test datapoints.
%               Of the form {index_train,index_test}.
%               If of the form [index], then just plots the subset index of the
%                       points in Features.
%               
% annotateParams, a structure containing all the information about how to annotate
%           data points. Fields can include:
%               - n, the number of data points to annotate
%               - userInput, 0: randomly selected datapoints, 1: user clicks to annotate datapoints
%               - fdim, 1x2 vector with width and height of time series as fraction of plot size
%               - maxL, maximum length of annotated time series
%               - textAnnotation: 'fileName', 'tsid', or 'none' to annotate this data
%               - cmap, a cell of colors, with elements for each group
%               - theMarkerSize, a custom marker size
%               - theLineWidth: line width for annotated time series
% 
% showDistr, if 1 (default), plots marginal density estimates for each variable
%                   (above and to the right of the plot), otherwise set to 0.
%                   
% classMethod, can select a classifier to fit to the different classes (e.g.,
%               'linclass' for a linear classifier).
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

function TS_plot_2d(Features,DataInfo,trainTest,annotateParams,showDistr,classMethod)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------

% Features should be a Nx2 vector of where to plot each of the N data objects in the two-dimensional space
if nargin < 1
    error('You must provide two-dimensional feature vectors for the data.')
end

% DataInfo should be a structure array with all the information about the data (same length as Features)
% Group should be a field in this structure array

if nargin < 3 || isempty(trainTest)
    trainTest = {};
end

if nargin < 4 || isempty(annotateParams)
    annotateParams = struct('n',0); % don't annotate
end

% By default, plot kernel density estimates above and on the side of the plot:
if nargin < 5 || isempty(showDistr)
    showDistr = 1;
end

if nargin < 6 || isempty(classMethod)
    classMethod = 'linclass';
end

makeFigure = 1; % default is to plot on a brand new figure('color','w')

% ------------------------------------------------------------------------------
%% Load data
% ------------------------------------------------------------------------------
% Data is not loaded, now it must be provided

labels = DataInfo.labels; % Feature labels
if isstruct(annotateParams) || annotateParams > 0
    DataLabels = DataInfo.DataLabels; % We need data labels
end
GroupNames = DataInfo.GroupNames;
GroupIndices = DataInfo.GroupIndices;
TimeSeriesData = DataInfo.TimeSeriesData;
numGroups = length(GroupNames);

% ------------------------------------------------------------------------------
%% Subset
% ------------------------------------------------------------------------------
% Only use a subset of the full matrix
if (length(trainTest)==1 || ~iscell(trainTest))
    if iscell(trainTest)
        rss = trainTest{1}; % row subset
    else
        rss = trainTest;
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
    trainTest = []; % make empty so don't plot trainTest groups later
end

% ------------------------------------------------------------------------------
% Compute classification rates
% ------------------------------------------------------------------------------
classRate = zeros(3,1); % classRate1, classRate2, classRateboth
groupLabels = BF_ToGroup(GroupIndices); % Convert GroupIndices to group form
switch classMethod
    case 'linclass'
        kfold = 10;
        numRepeats = 5;
        classify_fn = @(XTrain,yTrain,Xtest,ytest) ...
                        sum(ytest == classify(Xtest,XTrain,yTrain,'linear'))/length(ytest);
        try
            classRate(1) = mean(classify_fn(Features(:,1),groupLabels,Features(:,1),groupLabels));
            classRate(2) = mean(classify_fn(Features(:,2),groupLabels,Features(:,2),groupLabels));
            classRate(3) = mean(classify_fn(Features(:,1:2),groupLabels,Features(:,1:2),groupLabels));
        catch emsg
            fprintf(1,'%s\n',emsg.message);
        end
        fprintf(1,'Linear in-sample classification rates computed\n');
    case {'knn','knn_matlab'}
        k = 3;
        numRepeats = 10;
        try
            lf1 = TSQ_cfnerr('knn',k,Features(:,1),groupLabels,[],{'kfold',10,numRepeats});
            classRate(1,:) = [mean(lf1),std(lf1)];
            lf2 = TSQ_cfnerr('knn',k,Features(:,2),groupLabels,[],{'kfold',10,numRepeats});
            classRate(2,:) = [mean(lf2),std(lf2)];
            lf3 = TSQ_cfnerr('knn',k,Features(:,1:2),groupLabels,[],{'kfold',10,numRepeats});
            classRate(3,:) = [mean(lf3),std(lf3)];
        catch emsg
            fprintf(1,'%s\n',emsg);
        end
end


% ------------------------------------------------------------------------------
%% Plot
% ------------------------------------------------------------------------------
if makeFigure % can set extras.makeFigure = 0 to plot within a given setting
    f = figure('color','w'); box('on'); % white figure
    f.Position = [f.Position(1), f.Position(2), 600, 550];
end

% Set colors
if isstruct(annotateParams) && isfield(annotateParams,'cmap')
    if ischar(annotateParams.cmap)
        groupColors = BF_getcmap(annotateParams.cmap,numGroups,1);
    else
        groupColors = annotateParams.cmap; % specify the cell itself
    end
else
    if numGroups < 10
        groupColors = BF_getcmap('set1',numGroups,1);
    elseif numGroups <= 12
        groupColors = BF_getcmap('set3',numGroups,1);
    elseif numGroups <= 22
        groupColors = [BF_getcmap('set1',numGroups,1); ...
                    BF_getcmap('set3',numGroups,1)];
    elseif numGroups <= 50
        groupColors = mat2cell(jet(numGroups),ones(numGroups,1));
    else
        error('There aren''t enough colors in the rainbow to plot this many classes!')
    end
end
if (numGroups == 1)
    groupColors = {'k'}; % Just use black...
    % groupColors = 'rainbow'; % Just use black...
end


% ------------------------------------------------------------------------------
%% Plot distributions
% ------------------------------------------------------------------------------
if showDistr
    subplot(4,4,1:3); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:numGroups
        fr = plot_ks(Features(GroupIndices{i},1),groupColors{i},0);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'ylim',[minn,maxx]);
    
    subplot(4,4,[8,12,16]); hold on; box('on')
    maxx = 0; minn = 100;
    for i = 1:numGroups
        fr = plot_ks(Features(GroupIndices{i},2),groupColors{i},1);
        maxx = max([maxx,fr]); minn = min([minn,fr]);
    end
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    set(gca,'xlim',[minn,maxx]);
end


% ------------------------------------------------------------------------------
%% Set up a 2D plot
% ------------------------------------------------------------------------------
if showDistr
    subplot(4,4,[5:7,9:11,13:15]); box('on');
end
hold on;

if isfield(annotateParams,'theMarkerSize');
    theMarkerSize = annotateParams.theMarkerSize; % specify custom marker size
else
    if isempty(trainTest)
        theMarkerSize = 12; % Marker size for '.'
    else
        theMarkerSize = 5; % Marker size for 'o' and 's'
    end
end

if isempty(trainTest)
    for i = 1:numGroups
        plot(Features(GroupIndices{i},1),Features(GroupIndices{i},2),'.','color',groupColors{i},'MarkerSize',theMarkerSize)
    end
else % Plot training and test data differently
    for j = 1:length(trainTest)
        for i = 1:numGroups
            if (j==1)
                % Training data
                plot(Features(intersect(GroupIndices{i},trainTest{j}),1),Features(intersect(GroupIndices{i},trainTest{j}),2),...
                        'ok','MarkerFaceColor',groupColors{i},'MarkerSize',theMarkerSize)
            else
                % Test data
                plot(Features(intersect(GroupIndices{i},trainTest{j}),1),Features(intersect(GroupIndices{i},trainTest{j}),2),...
                        'sk','MarkerFaceColor',groupColors{i},'MarkerSize',theMarkerSize)
            end
        end
    end
end

% ------------------------------------------------------------------------------
%% Plot a classify boundary?
% ------------------------------------------------------------------------------
didClassify = 0;
if (numGroups == 2) && strcmp(classMethod,'linclass');
    modeorder = 'linear'; % or 'quadratic'
    
    xlim = get(gca,'XLim'); ylim = get(gca,'YLim');
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
    
    % Label that classification was performed
    didClassify = 1;
end


% ------------------------------------------------------------------------------
%% Label Axes
% ------------------------------------------------------------------------------
if didClassify
    title(sprintf('Combined classification rate (%s) = %.2f%%',classMethod, ...
                    round(classRate(3,1)*100)),'interpreter','none');
    labelText = cell(2,1);
    for i = 1:2
        labelText{i} = sprintf('%s (acc = %.2f %%)',labels{i}, ...
                                round(classRate(i,1)*100)); %,round(classRate(i,2)*100));
    end
else
    labelText = labels;
end

xlabel(labelText{1},'interpreter','none')
ylabel(labelText{2},'interpreter','none')

% Set Legend
if numGroups > 1
    if isempty(trainTest)
        legs = cell(numGroups,1);
        for i = 1:numGroups
            legs{i} = sprintf('%s (%u)',GroupNames{i},length(GroupIndices{i}));
        end
    else
        legs = cell(numGroups*2,1);
        for i = 1:numGroups
            legs{i} = sprintf('%s train (%u)',GroupNames{i},length(intersect(GroupIndices{i},trainTest{1})));
            legs{numGroups+i} = sprintf('%s test (%u)',GroupNames{i},length(intersect(GroupIndices{i},trainTest{2})));
        end
    end
    legend(legs);
end

% ------------------------------------------------------------------------------
%% Annotate time-series data
% ------------------------------------------------------------------------------
if isempty(TimeSeriesData)
    % Only attempt to annotate if time-series data is provided
    return
end
% Set parameters
if isfield(annotateParams,'maxL')
    maxL = annotateParams.maxL;
else
    maxL = 300; % length of annotated time series segments
end
if isfield(annotateParams,'userInput')
    userInput = annotateParams.userInput;
else
    userInput = 1; % user input points rather than randomly chosen
end
if isfield(annotateParams,'fdim')
    fdim = annotateParams.fdim;
else
    fdim = [0.30,0.08]; % width, height
end
if isfield(annotateParams,'textAnnotation')
    textAnnotation = annotateParams.textAnnotation; % 'fileName','tsid','none'
else
    textAnnotation = 'none'; % no annotations by default
end
if isfield(annotateParams,'theLineWidth')
    theLineWidth = annotateParams.theLineWidth;
else
    theLineWidth = 0.8; % % line width for time series
end
numAnnotations = annotateParams.n;

pxlim = get(gca,'xlim'); % plot limits
pylim = get(gca,'ylim'); % plot limits
pwidth = diff(pxlim); % plot width
pheight = diff(pylim); % plot height
alreadyPicked = zeros(numAnnotations,2); % record those already picked
plotCircle = 1; % magenta circle around annotated points
myColors = [BF_getcmap('set1',5,1);BF_getcmap('dark2',6,1)];

% Produce xy points
xy = cell(numGroups,1);
for i = 1:numGroups
    xy{i} = [Features(GroupIndices{i},1),Features(GroupIndices{i},2)];
end

% Don't use user input to select points to annotate: instead they are selected randomly
if ~userInput
    if numAnnotations == length(DataLabels) % annotate all
        fprintf(1,'Annotate all!\n')
        for j = 1:numAnnotations
            theGroup = find(cellfun(@(x)ismember(j,x),GroupIndices));
            alreadyPicked(j,1) = theGroup;
            alreadyPicked(j,2) = find(GroupIndices{theGroup}==j);
        end
    else
        alreadyPicked(:,1) = round(linspace(1,numGroups,numAnnotations));
        randPerms = cellfun(@(x)randperm(length(x)),GroupIndices,'UniformOutput',0);
        counters = ones(numGroups,1);
        for j = 1:numAnnotations
            alreadyPicked(j,2) = randPerms{alreadyPicked(j,1)}(counters(alreadyPicked(j,1))); % random element of the group
            counters(alreadyPicked(j,1)) = counters(alreadyPicked(j,1))+1;
        end
    end
end

% ------------------------------------------------------------------------------
% Go through and annotate selected points
fprintf(1,['Annotating time series segments to %u points in the plot, ' ...
                        'displaying %u samples from each...\n'],numAnnotations,maxL);
for j = 1:numAnnotations
    if userInput
        point = ginput(1);
        iplot = ClosestPoint_ginput(xy,point); % find closest actual point to input point
        theGroup = iplot(1); % want this group
        itsme = iplot(2); % and this index
        alreadyPicked(j,:) = [theGroup,itsme];
    else
        theGroup = alreadyPicked(j,1);
        itsme = alreadyPicked(j,2);
    end
    
    if (j > 1) && any(sum(abs(alreadyPicked(1:j-1,:) - repmat(alreadyPicked(j,:),j-1,1)),2)==0)
        % Same one has already been picked, don't plot it again
        continue
    end
    
    plotPoint = xy{theGroup}(itsme,:);
    theDataLabel = DataLabels{GroupIndices{theGroup}(itsme)}; % fileName of timeseries to plot
    timeSeriesSegment = TimeSeriesData{GroupIndices{theGroup}(itsme)}; % fileName of timeseries to plot
    if ~isempty(maxL)
        timeSeriesSegment = timeSeriesSegment(1:min(maxL,end));
    end
    
    % Plot a circle around the annotated point:
    if numGroups==1
        % cycle through rainvow colors sequentially:
        groupColors{1} = myColors{rem(j,length(myColors))};
    end
    if plotCircle
        plot(plotPoint(1),plotPoint(2),'o','MarkerEdgeColor',groupColors{theGroup},...
                            'MarkerFaceColor',brighten(groupColors{theGroup},0.5));
    end
    
    % Add text annotations:
    switch textAnnotation
    case 'fileName'
        % Annotate text with names of datapoints:
        text(plotPoint(1),plotPoint(2)-0.01*pheight,theDataLabel,...
                    'interpreter','none','FontSize',8,...
                    'color',brighten(groupColors{theGroup},-0.6));
    case 'tsid'
        % Annotate text with ts_id:
        text(plotPoint(1),plotPoint(2)-0.01*pheight,...
                num2str(ts_ids_keep(GroupIndices{theGroup}(itsme))),...
                    'interpreter','none','FontSize',8,...
                    'color',brighten(groupColors{theGroup},-0.6));
    end
    
    % Adjust if annotation goes off axis x-limits
    px = plotPoint(1)+[-fdim(1)*pwidth/2,+fdim(1)*pwidth/2];
    if px(1) < pxlim(1), px(1) = pxlim(1); end % can't plot off left side of plot
    if px(2) > pxlim(2), px(1) = pxlim(2)-fdim(1)*pwidth; end % can't plot off right side of plot
    
    % Adjust if annotation goes above maximum y-limits
    py = plotPoint(2)+[0,fdim(2)*pheight];
    if py(2) > pylim(2)
        py(1) = pylim(2)-fdim(2)*pheight;
    end
    
    % Annotate the time series
    plot(px(1)+linspace(0,fdim(1)*pwidth,length(timeSeriesSegment)),...
            py(1)+fdim(2)*pheight*(timeSeriesSegment-min(timeSeriesSegment))/(max(timeSeriesSegment)-min(timeSeriesSegment)),...
                '-','color',groupColors{theGroup},'LineWidth',theLineWidth);
end



% ------------------------------------------------------------------------------
function [fr, xr] = plot_ks(v,c,swap)
    % Vector v is the vector of a given group
    % c is the color
    [f, x] = ksdensity(v(~isnan(v)),linspace(min(v),max(v),1000),'function','pdf');
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
% ------------------------------------------------------------------------------
	
end