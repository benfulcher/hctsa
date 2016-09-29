function TS_plot_timeseries(whatData,numPerGroup,whatTimeSeries,maxLength,plotOptions)
% TS_plot_timeseries    Plots examples of time series in an hctsa analysis.
%
%---INPUTS:
% whatData, The hctsa data to load information from (cf. TS_LoadData)
%
% numPerGroup, If plotting groups, plots this many examples per group
%
% whatTimeSeries, Can provide indices to plot that subset, a keyword to plot
%                   matches to the keyword, 'all' to plot all, or an empty vector
%                   to plot default groups in TimeSeries.Group
%
% maxLength, the maximum number of samples of each time series to plot
%
% plotOptions, additional plotting options as a structure

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

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
% Get time-series data from where? ('norm' by default)
if nargin < 1 || isempty(whatData)
    whatData = 'norm';
end

if nargin < 2 || isempty(numPerGroup)
    % Default: plot 10 time series per group
    numPerGroup = 10;
end

% Can specify a reduced set of time series by keyword
if nargin < 3
    whatTimeSeries = '';
end

if nargin < 4
    % Maximum length of time series to display (otherwise crops)
    % If empty, displays all of all time series
    maxLength = [];
end

if nargin < 5
	plotOptions = [];
end

% ------------------------------------------------------------------------------
% Evaluate any custom plotting options specified in the structure plotOptions
% ------------------------------------------------------------------------------

if isstruct(plotOptions) && isfield(plotOptions,'displayTitles')
    displayTitles = plotOptions.displayTitles;
else
    % Show titles -- removing them allows more to be fit into plot
    displayTitles = 1; % show titles by default
end
if isstruct(plotOptions) && isfield(plotOptions,'howToFilter')
    howToFilter = plotOptions.howToFilter;
else
    howToFilter = 'evenly'; % by default
end
% Specify the colormap to use
if isstruct(plotOptions) && isfield(plotOptions,'colorMap')
    colorMap = plotOptions.colorMap;
else
    colorMap = ''; % choose automatically using GiveMeColors
end
% Specify whether to make a free-form plot
if isstruct(plotOptions) && isfield(plotOptions,'plotFreeForm')
    plotFreeForm = plotOptions.plotFreeForm;
else
    plotFreeForm = 1; % do a normal subplotted figure
end
% Specify line width for plotting
if isstruct(plotOptions) && isfield(plotOptions,'LineWidth')
    lw = plotOptions.LineWidth;
else
    lw = 1; % do a normal subplotted figure
end
% Determine whether to create a new figure
if isstruct(plotOptions) && isfield(plotOptions,'newFigure')
    newFigure = plotOptions.newFigure;
else
    newFigure = 1;
end

% ------------------------------------------------------------------------------
%% Load data
% ------------------------------------------------------------------------------
[~,TimeSeries] = TS_LoadData(whatData);

% ------------------------------------------------------------------------------
%% Get group indices:
% ------------------------------------------------------------------------------
if (isempty(whatTimeSeries) || strcmp(whatTimeSeries,'grouped')) && isfield(TimeSeries,'Group');
    % Use default groups
    groupIndices = BF_ToGroup([TimeSeries.Group]);
    fprintf(1,'Plotting from %u groups of time series from file.\n',length(groupIndices));
elseif isempty(whatTimeSeries) || strcmp(whatTimeSeries,'all')
    % Nothing specified but no groups assigned, or specified 'all': plot from all time series
    groupIndices = {1:length(TimeSeries)};
elseif ischar(whatTimeSeries)
    % Just plot the specified group
    % First load group names:
    if isstruct(whatData)
        groupNames = whatData.groupNames;
    else
        load(theFile,'groupNames');
    end
    a = strcmp(whatTimeSeries,groupNames);
    groupIndices = {find([TimeSeries.Group]==find(a))};
    fprintf(1,'Plotting %u time series matching group name ''%s''\n',length(groupIndices{1}),whatTimeSeries);
else % Provided a custom range as a vector
    groupIndices = {whatTimeSeries};
    fprintf(1,'Plotting the %u time series matching indices provided\n',length(whatTimeSeries));
end
numGroups = length(groupIndices);

% ------------------------------------------------------------------------------
%% Do the plotting
% ------------------------------------------------------------------------------
% Want numPerGroup from each time series group
iPlot = zeros(numGroups*numPerGroup,1);
classes = zeros(numGroups*numPerGroup,1);
nhere = zeros(numGroups,1);
groupSizes = cellfun(@length,groupIndices);

for i = 1:numGroups
    % filter down to numPerGroup if too many in group, otherwise plot all in
    % group
    switch howToFilter
        case 'firstcome'
            % just plot first in group (useful when ordered by closeness to
            % cluster centre)
            jj = (1:min(numPerGroup,groupSizes(i)));

        case 'evenly'
            % Plot evenly spaced through the given ordering
            jj = unique(round(linspace(1,groupSizes(i),numPerGroup)));

        case 'rand'
            % select ones to plot at random
            if groupSizes(i) > numPerGroup
                jj = randperm(groupSizes(i)); % randomly selected
                if length(jj) > numPerGroup
                    jj = jj(1:numPerGroup);
                end
            else
                jj = (1:min(numPerGroup,groupSizes(i))); % retain order if not subsampling
            end
    end
    nhere(i) = length(jj); % could be less than numPerGroup if a smaller group
    rh = sum(nhere(1:i-1))+1:sum(nhere(1:i)); % range here
    iPlot(rh) = groupIndices{i}(jj);
    classes(rh) = i;
end

% Summarize time series chosen to plot
rKeep = (iPlot > 0);
classes = classes(rKeep);
iPlot = iPlot(rKeep); % contains all the indicies of time series to plot (in order)
numToPlot = length(iPlot);

%-------------------------------------------------------------------------------
fprintf(1,'Plotting %u (/%u) time series from %u classes\n', ...
                    numToPlot,sum(cellfun(@length,groupIndices)),numGroups);

if isnumeric(colorMap)
    % Specified a custom colormap as a matrix
    theColors = mat2cell(colorMap);
elseif iscell(colorMap)
    % Specified a custom colormap as a cell of colors
    theColors = colorMap;
else
    % Get some colormap using GiveMeColors
    theColors = GiveMeColors(numGroups);
end

% ------------------------------------------------------------------------------
% Only create a new figure if required
if newFigure
    figure('color','w');
end

Ls = zeros(numToPlot,1); % length of each plotted time series
if plotFreeForm
    % FREEFORM: make all within a single plot with text labels
    ax = gca;
    ax.Box = 'on';
    hold(ax,'on');

	yr = linspace(1,0,numToPlot+1);
    inc = abs(yr(2)-yr(1)); % size of increment
    yr = yr(2:end);
	ls = zeros(numToPlot,1); % lengths of each time series
	if isempty(maxLength)
		for i = 1:numToPlot
			ls(i) = length(TimeSeries(iPlot(i)).Data);
		end
		maxN = max(ls); % maximum length of all time series to plot
	else
		maxN = maxLength;
	end

	for i = 1:numToPlot
	    fn = TimeSeries(iPlot(i)).Name; % the name of the time series
	    kw = TimeSeries(iPlot(i)).Keywords; % the keywords
	    x = TimeSeries(iPlot(i)).Data;
	    N0 = length(x);
		if ~isempty(maxN) && (N0 > maxN)
			% specified a maximum length of time series to plot
            sti = randi(N0-maxN,1);
			x = x(sti:sti+maxN-1); % subset random segment
            N = length(x);
        else
            N = N0; % length isn't changing
        end
		xx = (1:N) / maxN;
		xsc = yr(i) + 0.8*(x-min(x))/(max(x)-min(x)) * inc;

        if numGroups==1 && (length(theColors)==numToPlot)
            % Plot as per a set of colors provided:
            colorNow = theColors{i};
        else % plot by group color (or all black for 1 class)
            colorNow = theColors{classes(i)};
        end
        plot(xx,xsc,'-','color',colorNow,'LineWidth',lw)

        % Annotate text labels
		if displayTitles
			theTit = sprintf('{%u} %s [%s] (%u)',TimeSeries(iPlot(i)).ID,fn,kw,N0);
			text(0.01,yr(i)+0.9*inc,theTit,'interpreter','none','FontSize',8)
	    end
	end

    % Set up axes:
    ax.XTick = linspace(0,1,3);
    ax.XTickLabel = round(linspace(0,maxN,3));
    ax.YTick = [];
    ax.YTickLabel = {};
	ax.XLim = [0,1]; % Don't let the axes annoyingly slip out
    xlabel('Time (samples)')

else
    % i.e., NOT a FreeForm plot:
	for i = 1:numToPlot
	    subplot(numToPlot,1,i)
	    fn = TimeSeries(iPlot(i)).Name; % the filename
	    kw = TimeSeries(iPlot(i)).Keywords; % the keywords
	    x = TimeSeries(iPlot(i)).Data;
	    N = length(x);

	    % Prepare text for the title
		if displayTitles
			startBit = sprintf('{%u} %s [%s]',TimeSeries(iPlot(i)).ID,fn,kw);
	    end

	    % Plot the time series
	    if isempty(maxLength)
	        % no maximum length specified: plot the whole time series
	        plot(x,'-','color',theColors{classes(i)})
	        Ls(i) = N;
	        if displayTitles
	            title([startBit ' (' num2str(N) ')'],'interpreter','none','FontSize',8);
	        end
	    else
	        % Specified a maximum length of time series to plot: maxLength
	        if N <= maxLength
	            plot(x,'-','color',theColors{classes(i)});
	            Ls(i) = N;
	            if displayTitles
	                title([startBit ' (' num2str(N) ')'],'interpreter','none','FontSize',8);
	            end
	        else
	            sti = randi(N-maxLength,1);
	            plot(x(sti:sti+maxLength),'-','color',theColors{classes(i)}) % plot a random maxLength-length portion of the time series
	            Ls(i) = maxLength;
	            if displayTitles
	                title([startBit ' (' num2str(N) ' :: ' num2str(sti) '-' num2str(sti+maxLength) ')'],...
                                'interpreter','none','FontSize',8);
	            end
	        end
	    end
	    set(gca,'YTickLabel','');
	    if i~=numToPlot
	        set(gca,'XTickLabel','','FontSize',8) % put the ticks for the last time series
        else % label the axis
            xlabel('Time (samples)')
        end
	end

	% Set all xlims so that they have the same x-axis limits
	for i = 1:numToPlot
	    ax = subplot(numToPlot,1,i);
        ax.XLim = [1,max(Ls)];
	end
end


end
