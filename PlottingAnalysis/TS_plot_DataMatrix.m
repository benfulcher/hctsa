function TS_plot_DataMatrix(varargin)
% TS_plot_DataMatrix   Plot the data matrix.
%
%---EXAMPLE USAGE:
% TS_plot_DataMatrix; % plots clustered data if it exists
% TS_plot_DataMatrix('whatData','norm'); % plots normalized data
%
%---INPUTS:
% whatData: specify 'norm' for normalized data in HCTSA_N.mat, 'cl' for clustered
%         data in HCTSA_cl.mat (default), or specify a filename to load data
%         from that file.
% addTimeSeries: set to 1 to annotate time series segments to the left of the data matrix
% timeSeriesLength: length of time-series annotations (number of samples)
% colorGroups: set to 1 to color time-series groups with different colormaps
% customColorMap: use a custom color map (name to match an option in BF_getcmap)
% colorNaNs: whether to plot rectangles over special-values in the matrix (default: 1)
% customOrder: reorder rows and columns according to provided permutation vectors
%
%---OUTPUT:
% Produces a colormap plot of the data matrix with time series as rows and
%   operations as columns.

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
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
inputP = inputParser;

% whatDataFile
default_whatData = 'cl';
check_whatData = @(x) ischar(x) || isstruct(x);
addOptional(inputP,'whatData',default_whatData,check_whatData);

% addTimeSeries, annotates time series segments to the side of the plot
default_addTimeSeries = 1;
check_addTimeSeries = @(x) isnumeric(x) && (x==0 || x==1);
addOptional(inputP,'addTimeSeries',default_addTimeSeries,check_addTimeSeries);

% timeSeriesLength, length of time-series annotations to the left of the main plot
default_timeSeriesLength = 100;
addOptional(inputP,'timeSeriesLength',default_timeSeriesLength,@isnumeric);

% colorGroups, color groups of time series differently:
default_colorGroups = 0;
check_colorGroups = @(x) (x==0 || x==1);
addOptional(inputP,'colorGroups',default_colorGroups,check_colorGroups);

% groupReorder, reorder within groups of time series:
default_groupReorder = 0;
check_groupReorder = @(x) (x==0 || x==1);
addOptional(inputP,'groupReorder',default_groupReorder,check_groupReorder);

% custom color map, customColorMap
default_customColorMap = 'redyellowblue';
addOptional(inputP,'customColorMap',default_customColorMap,@ischar);

% colorNaNs
default_colorNaNs = 1;
check_colorNaNs = @(x) (x==0 || x==1);
addOptional(inputP,'colorNaNs',default_colorNaNs,check_colorNaNs);

% customOrder (reorder before plotting)
default_customOrder = {[],[]};
check_customOrder = @(x)iscell(x) && length(x)==2;
addOptional(inputP,'customOrder',default_customOrder,check_customOrder);

%% Parse inputs:
parse(inputP,varargin{:});

% Make variables from results of input parser:
addTimeSeries = inputP.Results.addTimeSeries;
whatData = inputP.Results.whatData;
timeSeriesLength = inputP.Results.timeSeriesLength;
colorGroups = inputP.Results.colorGroups;
groupReorder = inputP.Results.groupReorder;
customColorMap = inputP.Results.customColorMap;
colorNaNs = inputP.Results.colorNaNs;
customOrder = inputP.Results.customOrder;
clear inputP;

% --------------------------------------------------------------------------
%% Read in the data
% --------------------------------------------------------------------------
% You always want to retrieve and plot the clustered data if it exists
getClustered = 1;
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData,getClustered);

[numTS, numOps] = size(TS_DataMat); % size of the data matrix

% ------------------------------------------------------------------------------
%% Reorder according to customOrder
% ------------------------------------------------------------------------------
if ~isempty(customOrder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified.\n');
	TS_DataMat = TS_DataMat(customOrder{1},:);
    TimeSeries = TimeSeries(customOrder{1});
end

if ~isempty(customOrder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified.\n');
	TS_DataMat = TS_DataMat(:,customOrder{2});
    Operations = Operations(customOrder{2});
end

%-------------------------------------------------------------------------------
% Check group information
%-------------------------------------------------------------------------------
if isfield(TimeSeries,'Group')
	timeSeriesGroups = [TimeSeries.Group];
	numClasses = max(timeSeriesGroups);
else
	timeSeriesGroups = [];
end
if colorGroups==1
	if ~isempty(timeSeriesGroups)
	    fprintf(1,'Coloring groups of time series...\n');
	else
	    warning('No group information found')
	    colorGroups = 0;
	end
end

%-------------------------------------------------------------------------------
% Reorder according to groups
%-------------------------------------------------------------------------------
if groupReorder
	if isempty(timeSeriesGroups)
		warning('Cannot reorder by time series group; no group information found')
	else
	    [~,ixData] = sort(timeSeriesGroups,'ascend');
	    dataMatReOrd = TS_DataMat(ixData,:);
	    ixAgain = ixData;
	    for i = 1:numClasses
	        isGroup = [TimeSeries(ixData).Group]==i;
	        ordering = BF_ClusterReorder(dataMatReOrd(isGroup,:),'euclidean','average');
	        istmp = ixData(isGroup);
	        ixAgain(isGroup) = istmp(ordering);
	    end
	    ixData = ixAgain; % set ordering to ordering within groups
	    TimeSeries = TimeSeries(ixData);
	    TS_DataMat = TS_DataMat(ixData,:);
	    timeSeriesGroups = timeSeriesGroups(ixData);
	end
end

% --------------------------------------------------------------------------
%% Prepare data matrix for plotting
% --------------------------------------------------------------------------
if colorGroups
    gi = BF_ToGroup(timeSeriesGroups);

    numGroups = length(gi);

    % Add a group for unlabelled data items if they exist
    if sum(cellfun(@length,gi)) < numTS
        % Add an unlabelled class
        gi0 = gi;
        gi = cell(numGroups+1,1);
        for i = 1:numGroups
            gi{i} = gi0{i};
        end
        clear gi0;
        gi{end} = setxor(1:numTS,vertcat(gi{:}));
        numGroups = numGroups + 1;
    end

    fprintf(1,'Coloring data according to %u groups\n',numGroups);

    % Change range of TS_DataMat to make use of new colormap appropriately
    ff = 0.9999999;
    squashMe = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    TS_DataMat = squashMe(TS_DataMat);
    for jo = 1:numGroups
        TS_DataMat(gi{jo},:) = squashMe(TS_DataMat(gi{jo},:)) + jo - 1;
    end
else
    numGroups = 0;
end

% --------------------------------------------------------------------------
%% Set the colormap
% --------------------------------------------------------------------------
if numGroups <= 1
    numColorMapGrads = 6; % number of gradations in each set of colourmap
    if strcmp(customColorMap,'redyellowblue');
        customColorMap = flipud(BF_getcmap('redyellowblue',numColorMapGrads,0));
    else
        customColorMap = gray(numColorMapGrads);
    end
elseif numGroups ==2
	% Special case to make a nice red and blue one
	customColorMap = [flipud(BF_getcmap('blues',9,0));flipud(BF_getcmap('reds',9,0))];
else
	% Use the same colors as GiveMeColors, but add brightness gradations to indicate magnitude
    numColorMapGrads = 20; % number of brightness gradations in each set of colourmap
    colormapBase = GiveMeColors(numGroups);
    customColorMap = [];
    for i = 1:numGroups
        customColorMap = [customColorMap; BF_MakeBrightenedColorMap(colormapBase{i},numColorMapGrads)];
    end
end

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
f = figure('color','w');

if addTimeSeries
    %  First make an additional plot to include time-series subsegments
    ax1 = subplot(1,5,1);
    hold(ax1,'on');
    ax1.Box = 'on';
    ax1.YTick = (1:numTS);
    ax1.YTickLabel = {TimeSeries.Name};
    ax1.YLim = [0.5,numTS+0.5];
    ax1.XLim = [1,timeSeriesLength];
    xlabel('Time (samples)');
    ax1.TickLabelInterpreter = 'none';
    NormMinMax = @(x) (x-min(x))/(max(x)-min(x));
    for j = 1:numTS
        % Plot a segment from each time series, up to a maximum length of
        % timeSeriesLength samples (which is set as an input to the function)
        tsData = TimeSeries(j).Data;
        lengthHere = min(timeSeriesLength,length(tsData));
        tsData = tsData(1:lengthHere);
        plot(1:lengthHere,j-0.5+NormMinMax(tsData),'-k');
        if j < numTS
            plot([1,timeSeriesLength],(j+0.5)*ones(2,1),':k')
        end
    end

    ax2 = subplot(1,5,2:5);
    hold(ax2,'on');
    ax.Box = 'on';

    linkaxes([ax1,ax2],'y');
end

% ------------------------------------------------------------------------------
%% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
% (alternative is to use imagesc)
colormap(customColorMap)
imagesc(TS_DataMat);

% ------------------------------------------------------------------------------
% Superimpose colored rectangles over NaN values
% ------------------------------------------------------------------------------
if colorNaNs && any(isnan(TS_DataMat(:)))
    [theNaNs_i,theNaNs_j] = find(isnan(TS_DataMat));
    fprintf(1,['Superimposing black rectangles over all %u NaNs in ' ...
                                'the data matrix\n'],length(theNaNs_i));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end
end

% --------------------------------------------------------------------------
%% Format the axes
% --------------------------------------------------------------------------
% Axis labels:
ax2 = gca;
ax2.FontSize = 8; % small font size (often large datasets)
ax2.TickLabelInterpreter = 'none'; % Stop from displaying underscores as subscripts

% Rows: time series
ax2.YTick = 1:numTS;
ax2.YLim = [0.5,numTS+0.5];

% Columns: operations:
xlabel('Operations')
ax2.XLim = [0.5,numOps+0.5];
if numOps < 1000 % if too many operations, it's too much to list them all...
    ax2.XTick = 1:numOps;
    ax2.XTickLabel = {Operations.Name};
    ax2.XTickLabelRotation = 90;
end

% Add a color bar:
cB = colorbar('eastoutside');
cB.Label.String = 'Output';
if numGroups > 0
	cB.Ticks = 0.5:1:numGroups;
	cB.TickLabels = TS_GetFromData(whatData,'groupNames');
end

title(sprintf('Data matrix (%u x %u)',numTS,numOps))

if addTimeSeries
    ax2.YTickLabel = {};
    % Reposition tight
    ax2.Position = [0.4,0.1,0.55,0.85];
    ax1.Position = [0.05,0.1,0.2,0.85];
    % Tight on left:
    ax1.Position(1) = ax2.Position(1) - ax1.Position(3) - 0.01;
    % Both have the same y and height:
    ax1.Position(2) = ax2.Position(2); % y
    ax1.Position(4) = ax2.Position(4); % height
else
    % Rows -- time series need labels
    ylabel('Time series')
    ax2.YTickLabel = {TimeSeries.Name};
end

end
