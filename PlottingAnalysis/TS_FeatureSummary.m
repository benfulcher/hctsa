function TS_FeatureSummary(opID, whatData, annotateParams)
% TS_FeatureSummary   How a given feature behaves across a time-series dataset
%
% Plots the distribution of outputs of an operation across the given dataset
% and allows the user to annotate time series onto the plot to visualize
% how the operation is behaving.
%
%---INPUTS:
% opID, the operation ID to plot
% whatData, the data to visualize (HCTSA_loc.mat by default; cf. TS_LoadData)
% annotateParams, a structure of custom plotting options

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

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------

if nargin < 2 || isempty(whatData)
   whatData = 'loc'; % Visualize at unnormalized outputs
end

if nargin < 3 || isempty(annotateParams) % annotation parameters
    annotateParams = struct(); % take 10 annotatations
end

%-------------------------------------------------------------------------------
% Load in data:
%-------------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations,whatDataFile] = TS_LoadData(whatData);
numTimeSeries = length(TimeSeries);
theOp = ([Operations.ID]==opID);
dataVector = TS_DataMat(:,theOp); % the outputs of interest
notNaN = find(~isnan(dataVector));
dataVector = dataVector(notNaN); % remove bad values
theOperation = Operations(theOp);

if isempty(dataVector)
    error('No data for %s',Operations(theOp).Name);
end

%-------------------------------------------------------------------------------
% Sort out any custom plotting settings in the annotateParams structure
%-------------------------------------------------------------------------------
if isfield(annotateParams,'n')
    numAnnotate = annotateParams.n;
else
    numAnnotate = min(6,numTimeSeries);
end
if isfield(annotateParams,'maxL')
    maxL = annotateParams.maxL;
else
    maxL = 300;
end
if isfield(annotateParams,'fdim')
    fdim = annotateParams.fdim;
else
    fdim = [0.30,0.08]; % width, height
end
if isfield(annotateParams,'uinput')
    uinput = annotateParams.uinput;
else
    uinput = 1; % user clicks to annotate rather than randomly chosen
end
if isfield(annotateParams,'textann') % what text annotation to use
    textann = annotateParams.textann;
else
    textann = 'filename'; %'filename','tsid','none','length'
end
if isfield(annotateParams,'whereann') % what text annotation to use
    textann = annotateParams.whereann;
else
    whereann = 'onplot'; %'onplot', 'newplot'
end

%-------------------------------------------------------------------------------
%% Plot the kernel smoothed density
%-------------------------------------------------------------------------------
numPoints = 1000;
fig = figure('color','w'); box('on'); hold on
if strcmp(whereann,'newplot')
    subplot(3,1,[1,2]); box('on'); hold on
end
[f,x] = ksdensity(dataVector,linspace(min(dataVector),max(dataVector),numPoints),'function','pdf');
getIndex = @(m) find(x>=dataVector(m),1,'first');
ind = arrayfun(getIndex,1:length(dataVector));
plot(x,f,'k'); % the curve
plot(x(ind),f(ind),'.k','MarkerSize',8); % individual TimeSeries as points


%-------------------------------------------------------------------------------
%% Annotate time series:
%-------------------------------------------------------------------------------
% Vector of points in the space, xy:
xy = [x(ind)',f(ind)'];
% Want on a common scale for finding neighbors:
xy_std = std(xy);
xy_mean = mean(xy);
xy_zscore = zscore(xy);

pWidth = diff(get(gca,'xlim')); % plot width
pHeight = diff(get(gca,'ylim')); % plot height
alreadyPicked = zeros(numAnnotate,1); % record those already picked
lineWidth = 0.7; % line width for time series
pxlim = get(gca,'xlim'); % plot limits
pylim = get(gca,'ylim'); % plot limits
cblue = BF_getcmap('set1',2,1); cblue = cblue{2};

plotCircle = 1; % magenta circle around annotated points
if ~uinput % points to annotate are randomly picked
    rp = randperm(numTimeSeries);
    alreadyPicked = rp(1:numAnnotate);
end

xlabel(theOperation.Name,'Interpreter','none');
ylabel('Probability Density')
for j = 1:numAnnotate
    title(sprintf('%u points remaining to annotate',numAnnotate-j+1));

    if uinput % user input
        point = ginput(1);
        point_z = (point-xy_mean)./xy_std;
        iplot = BF_ClosestPoint_ginput(xy_zscore,point_z); % find closest actual point to input point
        iplot = iplot(2);
        alreadyPicked(j) = iplot;
    else
        iplot = alreadyPicked(j);
    end

    if j > 1 && any(alreadyPicked(1:j-1)==iplot) % same already been picked
        beep; continue; % don't plot this again
    end

    plotpoint = xy(iplot,:);
    iFF = notNaN(iplot); % index of dataVector
    fn = TimeSeries(iFF).Name; % filename of timeseries to plot
    ts = TimeSeries(iFF).Data; % data of time series to plot
    if strcmp(whereann,'onplot')
        % Add time series traces as little annotations on the distribution plot
        if plotCircle
            plot(plotpoint(1),plotpoint(2),'o','color',cblue); % plot magenta circle around target point
        end
        if strcmp(textann,'filename')
            % annotate text with filename:
            text(plotpoint(1),plotpoint(2)-0.01*pHeight,fn,'interpreter','none','FontSize',8);
        elseif strcmp(textann,'tsid');
            % annotate text with ts_id:
            text(plotpoint(1),plotpoint(2)-0.01*pHeight,num2str(ts_ids_keep(iplot)),'interpreter','none','FontSize',8);
        elseif strcmp(textann,'length')
            text(plotpoint(1),plotpoint(2)-0.01*pHeight,num2str(length(ts)),'interpreter','none','FontSize',8);
        end
        if ~isempty(maxL)
            ts = ts(1:min(maxL,end)); % crop the time series
        end

        % adjust if annotation goes off axis x-limits
        px = plotpoint(1)+[-fdim(1)*pWidth/2,+fdim(1)*pWidth/2];
        if px(1)<pxlim(1), px(1) = pxlim(1); end % can't plot off left side of plot
        if px(2)>pxlim(2), px(1) = pxlim(2)-fdim(1)*pWidth; end % can't plot off right side of plot

        % adjust if annotation goes above maximum y-limits
        py = plotpoint(2)+[0,fdim(2)*pHeight];
        if py(2)>pylim(2), py(1) = pylim(2)-fdim(2)*pHeight; end

        plot(px(1)+linspace(0,fdim(1)*pWidth,length(ts)),...
            py(1)+fdim(2)*pHeight*(ts-min(ts))/(max(ts)-min(ts)),...
            '-','color','k','LineWidth',lineWidth);

    else
        % Plot time series traces in a subplot below the distribution -- to
        % visualize one at a time
        subplot(3,1,[1,2]); hold on
        if plotCircle
            plot(plotpoint(1),plotpoint(2),'om'); % plot magenta circle around target point
        end
        subplot(3,1,3);
        if ~isempty(maxL)
            ts = ts(1:min(maxL,end)); % crop the time series
        end
        plot(BF_zscore(ts),'.-k');
        xlabel(sprintf('[%u] %s (%u)',TimeSeries(iplot).ID,TimeSeries(iplot).Name,length(ts)))

        subplot(3,1,[1,2])
        title(num2str(numAnnotate-j));
    end
end
title(theOperation.Name,'Interpreter','none');
xlabel('Outputs','Interpreter','none');

end
