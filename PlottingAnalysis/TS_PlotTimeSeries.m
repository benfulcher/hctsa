function TS_PlotTimeSeries(whatData,numPerGroup,whatTimeSeries,maxLength,plotOptions)
% TS_PlotTimeSeries    Plots examples of time series in an hctsa analysis.
%
%---INPUTS:
% whatData, The hctsa data to load information from (cf. TS_LoadData) [or can
%               specify a TimeSeries table directly]
%
% numPerGroup, If plotting groups, plots this many examples per group
%
% whatTimeSeries, Can provide indices to plot that subset (rows of the time series
%                   table, after any sorting by plotOptions.sortByLength), the name
%                   of a group in TimeSeries.Group to plot that group, a keyword to
%                   plot the time series matching the keyword, 'all' to plot all,
%                   or an empty vector (or 'grouped') to plot group information
%                   assigned to TimeSeries.Group
%
% maxLength, the maximum number of samples of each time series to plot
%
% plotOptions, additional plotting options as a structure, with optional fields:
%       displayTitles, whether to label each time series (default: true)
%       howToFilter, how to choose numPerGroup time series from a group: 'evenly'
%                   (default; evenly spaced through the ordering), 'firstcome' (the
%                   first in order), or 'rand' (a random set)
%       colorMap, a cell of colors, a matrix of colors (one row per group), or the
%                   name of a color map (cf. BF_GetColorMap) (default: GiveMeColors)
%       plotFreeForm, true (default): all time series in a single set of axes;
%                   false: one subplot per time series
%       LineWidth, line width of free-form plots (default: 1)
%       newFigure, whether to plot in a new figure (default: true)
%       sortByLength, whether to sort time series by length first (default: false)
%       carpetPlot, whether to plot a grayscale image of the time series instead
%                   (default: false)

% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
% Get time-series data from where? ('norm' by default)
if nargin < 1 || isempty(whatData)
    whatData = 'raw';
end

if nargin < 2 || isempty(numPerGroup)
    % Default: plot 10 time series per group
    numPerGroup = 10;
end
if numPerGroup <= 0
    error('numPerGroup must be positive');
end

% Can specify a reduced set of time series by keyword (or indices)
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

%-------------------------------------------------------------------------------
% Evaluate any custom plotting options specified in the structure plotOptions
%-------------------------------------------------------------------------------
if isstruct(plotOptions) && isfield(plotOptions,'displayTitles')
    displayTitles = plotOptions.displayTitles;
else
    % Show titles -- removing them allows more to be fit into plot
    displayTitles = true; % show titles by default
end
if isstruct(plotOptions) && isfield(plotOptions,'howToFilter')
    howToFilter = plotOptions.howToFilter;
else
    howToFilter = 'evenly'; % by default
    % 'firstcome' (first time series in order)
    % 'evenly' (evenly spaced across the ordering)
    % 'rand' (random set; picks different time series each time)
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
    plotFreeForm = true; % plot all time series in a single set of axes
end
% Specify line width for plotting
if isstruct(plotOptions) && isfield(plotOptions,'LineWidth')
    lw = plotOptions.LineWidth;
else
    lw = 1; % line width
end
% Determine whether to create a new figure
if isstruct(plotOptions) && isfield(plotOptions,'newFigure')
    newFigure = plotOptions.newFigure;
else
    newFigure = true;
end
% Sorting time series by length can help visualization
if isstruct(plotOptions) && isfield(plotOptions,'sortByLength')
    sortByLength = plotOptions.sortByLength;
else
    sortByLength = false;
end
% Carpet plot instead...?
if isstruct(plotOptions) && isfield(plotOptions,'carpetPlot')
    carpetPlot = plotOptions.carpetPlot;
else
    carpetPlot = false;
end

% ------------------------------------------------------------------------------
%% Load TimeSeries table
% ------------------------------------------------------------------------------
if istable(whatData)
    TimeSeries = whatData;
else
    TimeSeries = TS_GetFromData(whatData,'TimeSeries');
end
if ismember('Group',TimeSeries.Properties.VariableNames)
    TimeSeries.Group = removecats(TimeSeries.Group); % only count groups with members
end
if sortByLength
    TimeSeries = sortrows(TimeSeries,'Length','descend');
end

% ------------------------------------------------------------------------------
%% Get group indices:
% ------------------------------------------------------------------------------
if (isempty(whatTimeSeries) || strcmp(whatTimeSeries,'grouped')) && ismember('Group',TimeSeries.Properties.VariableNames)
    % Use default groups assigned by TS_LabelGroups
    groupIndices = BF_ToGroup(TimeSeries.Group);
    classLabels = categories(TimeSeries.Group);
    numGroups = length(classLabels);
    fprintf(1,'Plotting from %u groups of time series from file.\n',numGroups);
elseif isempty(whatTimeSeries) || (ischar(whatTimeSeries) && ismember(whatTimeSeries,{'all','grouped'}))
    % Nothing specified but no groups assigned, or specified 'all': plot from all time series
    groupIndices = {1:height(TimeSeries)};
    classLabels = {};
elseif ischar(whatTimeSeries)
    if ismember('Group',TimeSeries.Properties.VariableNames) && ...
                ismember(whatTimeSeries,categories(TimeSeries.Group))
        % Just plot the specified group
        groupIndices = {find(TimeSeries.Group==whatTimeSeries)};
        classLabels = {whatTimeSeries};
        fprintf(1,'Plotting %u time series matching group name ''%s''\n',length(groupIndices{1}),whatTimeSeries);
    else
        % Plot the time series with a matching keyword
        Keywords = SUB_cell2cellcell(TimeSeries.Keywords);
        isMatch = cellfun(@(x) ismember(whatTimeSeries,x),Keywords);
        if ~any(isMatch)
            error('No group or keyword named ''%s'' found',whatTimeSeries);
        end
        groupIndices = {find(isMatch)};
        classLabels = {};
        fprintf(1,'Plotting %u time series matching keyword ''%s''\n',length(groupIndices{1}),whatTimeSeries);
    end
else % Provided a custom range as a vector
    groupIndices = {whatTimeSeries};
    classLabels = {};
    fprintf(1,'Plotting the %u time series matching indices provided\n',length(whatTimeSeries));
end
numGroups = length(groupIndices);

%-------------------------------------------------------------------------------
%% Do the plotting
%-------------------------------------------------------------------------------
% Want to plot numPerGroup from each time-series group
iPlot = zeros(numGroups*numPerGroup,1);
classes = zeros(numGroups*numPerGroup,1);
nhere = zeros(numGroups,1);
groupSizes = cellfun(@length,groupIndices);

for i = 1:numGroups
    % filter down to numPerGroup if too many in group, otherwise plot all in group
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
            jj = sort(jj,'ascend'); % Order time series in ascending order of index
        otherwise
            error('Unknown filtering option ''%s''',howToFilter);
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
%-------------------------------------------------------------------------------

% Create a new figure if required
if newFigure
    figure('color','w');
end

% Set default for max length if unspecified (as max length of time series)
if isempty(maxLength)
    ls = TimeSeries.Length(iPlot);
    maxN = max(ls); % maximum length of all time series to plot
else
    maxN = maxLength;
end

% Carpet plot is just a grayscale visualization of many time series
if carpetPlot
    % Assemble time-series data matrix:
    X = nan(numToPlot,maxN);
    for i = 1:numToPlot
        x = TimeSeries.Data{iPlot(i)};
        L = min(maxN,length(x));
        X(i,1:L) = x(1:L);
        % NB: other plotting uses random subsegment when longer; here just plot first maxN samples
    end
    % Normalize:
    Xnorm = BF_NormalizeMatrix(X','maxmin')'; % linearly normalize across rows
    imagesc(Xnorm);

    colormap(gray)

    % Add filenames to axes:
    if displayTitles
        ax = gca;
        fn = TimeSeries.Name(iPlot); % the names of the time series
        ax.YTick = 1:numToPlot;
        ax.YTickLabel = fn;
        ax.TickLabelInterpreter = 'none';
    end
    xlabel('Time (samples)')
    return
end

%-------------------------------------------------------------------------------
% Set colormap
%-------------------------------------------------------------------------------
if isnumeric(colorMap)
    % Specified a custom colormap as a matrix
    theColors = num2cell(colorMap,2); % one RGB row per color
elseif iscell(colorMap)
    % Specified a custom colormap as a cell of colors
    theColors = colorMap;
elseif ischar(colorMap) && ~isempty(colorMap)
    % Specified the name of a colormap
    theColors = BF_GetColorMap(colorMap,numGroups,true);
else
    % Get some colormap using GiveMeColors
    theColors = GiveMeColors(numGroups);
end
if numGroups > 1 && length(theColors) < numGroups
    error('colorMap specifies %u colors but there are %u groups to plot',length(theColors),numGroups);
end

%-------------------------------------------------------------------------------
% Plot as conventional time series
%-------------------------------------------------------------------------------
Ls = zeros(numToPlot,1); % length of each plotted time series
if plotFreeForm
    % FREEFORM: make all within a single plot with text labels
    ax = gca;
    ax.Box = 'on';
    hold(ax,'on');

    yr = linspace(1,0,numToPlot+1);
    inc = abs(yr(2)-yr(1)); % size of increment
    yr = yr(2:end);

    pHandles = zeros(numToPlot,1); % keep plot handles
    for i = 1:numToPlot
        x = TimeSeries.Data{iPlot(i)};
        N0 = length(x);
        if ~isempty(maxN) && (N0 > maxN)
            % Specified a maximum length of time series to plot
            sti = 1; % randi(N0-maxN,1);
            x = x(sti:sti+maxN-1); % subset random segment
            N = length(x);
        else
            N = N0; % length isn't changing
        end
        xx = (1:N) / maxN;
        xRange = max(x) - min(x);
        if xRange==0
            xRange = 1; % a constant time series is plotted flat
        end
        xsc = yr(i) + 0.8*(x-min(x))/xRange * inc;

        if numGroups==1 && (length(theColors)==numToPlot)
            % Plot as per a set of colors provided:
            colorNow = theColors{i};
        else % plot by group color (or all black for 1 class)
            colorNow = theColors{classes(i)};
        end
        pHandles(i) = plot(xx,xsc,'-','color',colorNow,'LineWidth',lw);

        % Annotate text labels
        if displayTitles
            theTit = sprintf('{%u} %s [%s] (%u)',TimeSeries.ID(iPlot(i)),...
                    TimeSeries.Name{iPlot(i)},TimeSeries.Keywords{iPlot(i)},N0);
            text(0.01,yr(i)+0.9*inc,theTit,'interpreter','none','FontSize',8)
        end
    end

    % Legend:
    if ~isempty(classLabels)
        [~,b] = unique(classes);
        legend(pHandles(b),classLabels,'interpreter','none');
    end

    % Set up axes:
    ax.XTick = linspace(0,1,3);
    ax.XTickLabel = round(linspace(0,maxN,3));
    ax.YTick = [];
    ax.YTickLabel = {};
	ax.XLim = [0,1]; % Don't let the axes annoyingly slip out
    xlabel('Time (samples)')

else
    % NOT a 'free-form' plot:
    axHandles = gobjects(numToPlot,1);
    for i = 1:numToPlot
        ax = subplot(numToPlot,1,i);
        axHandles(i) = ax;
        fn = TimeSeries.Name{iPlot(i)}; % the filename
        kw = TimeSeries.Keywords{iPlot(i)}; % the keywords
        x = TimeSeries.Data{iPlot(i)};
        N = length(x);

        % Prepare text for the title
        if displayTitles
            startBit = sprintf('{%u} %s [%s]',TimeSeries.ID(iPlot(i)),fn,kw);
        end

        % Plot the time series
        if isempty(maxLength)
            % no maximum length specified: plot the whole time series
            plot(x,'-','color',theColors{classes(i)})
            Ls(i) = N;
            if displayTitles
                title(sprintf('%s (%u)',startBit,N),'interpreter','none','FontSize',8);
            end
        else
            % Specified a maximum length of time series to plot: maxLength
            if N <= maxLength
                plot(x,'-','color',theColors{classes(i)});
                Ls(i) = N;
                if displayTitles
                    title(sprintf('%s (%u)',startBit,N),'interpreter','none','FontSize',8);
                end
            else
                sti = randi(N-maxLength,1);
                plot(x(sti:sti+maxLength-1),'-','color',theColors{classes(i)}) % plot a random maxLength-length portion of the time series
                Ls(i) = maxLength;
                if displayTitles
                    title(sprintf('%s (%u :: %u-%u)',startBit,N,sti,sti+maxLength-1),...
                                'interpreter','none','FontSize',8);
                end
            end
        end
        ax.YTickLabel = '';
        if i~=numToPlot
            ax.XTickLabel = '';
            ax.FontSize = 8; % put the ticks for the last time series
        else % label the axis
            ax.XLabel.String = 'Time (samples)';
        end
    end

    % Set all subplots to have the same x-axis limits
    for i = 1:numToPlot
        axHandles(i).XLim = [1,max(Ls)];
    end
end

end
