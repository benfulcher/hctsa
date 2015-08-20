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
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>

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

% ------------------------------------------------------------------------------
%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
inputP = inputParser;

% whatDataFile
default_whatData = 'cl';
check_whatData = @(x)1;
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
customColorMap = inputP.Results.customColorMap;
colorNaNs = inputP.Results.colorNaNs;
customOrder = inputP.Results.customOrder;
clear inputP;

% --------------------------------------------------------------------------
%% Read in the data
% --------------------------------------------------------------------------
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData);

if colorGroups
    if isfield(TimeSeries,'Group')
        timeSeriesGroups = [TimeSeries.Group];
        fprintf(1,'Coloring groups of time series...\n');
    else
        warning('No group information found')
        colorGroups = 0;
    end
end

TimeSeries; % Just extract time series names
Operations; % Just extract operation names

[numTS, numOps] = size(TS_DataMat); % size of the data matrix

% ------------------------------------------------------------------------------
%% Reorder according to customOrder
% ------------------------------------------------------------------------------
if ~isempty(customOrder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified.\n');
	TS_DataMat = TS_DataMat(customOrder{1},:);
    TimeSeries = timeSeries(customOrder{1});
end

if ~isempty(customOrder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified.\n');
	TS_DataMat = TS_DataMat(:,customOrder{2});
    Operations = Operations(customOrder{2});
end

% --------------------------------------------------------------------------
%% Prepare data matrix for plotting
% --------------------------------------------------------------------------
numColorMapGrads = 6; % number of gradations in each set of colourmap

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
    if strcmp(customColorMap,'redyellowblue');
        customColorMap = flipud(BF_getcmap('redyellowblue',numColorMapGrads,0));
    else
        customColorMap = gray(numColorMapGrads);
    end
else
    customColorMap = colormap(BF_getcmap('blues',numColorMapGrads,0,1));
    if numGroups >= 2
        customColorMap = [customColorMap; BF_getcmap('greens',numColorMapGrads,0,1)];
    end
    if numGroups >= 3
        customColorMap = [customColorMap; BF_getcmap('oranges',numColorMapGrads,0,1)];
    end
    if numGroups >= 4
        customColorMap = [customColorMap; BF_getcmap('purples',numColorMapGrads,0,1)];
    end
    if numGroups >= 5
        customColorMap = [customColorMap; BF_getcmap('reds',numColorMapGrads,0,1)];
    end
    if numGroups >= 6
        customColorMap = [customColorMap; pink(numColorMapGrads)];
    end
    if numGroups >= 7
        customColorMap = [customColorMap; gray(numColorMapGrads)];
    end
    if numGroups >= 8
        customColorMap = [customColorMap; BF_getcustomColorMap('yelloworangered',numColorMapGrads,0,1)];
    end
    if numGroups >= 9
        customColorMap = [customColorMap; BF_getcmap('purplebluegreen',numColorMapGrads,0,1)];
    end
    if numGroups >= 10
        customColorMap = [customColorMap; BF_getcmap('yellowgreenblue',numColorMapGrads,0,1)];
    end
    if numGroups >= 11
        customColorMap = [customColorMap; BF_getcmap('purpleblue',numColorMapGrads,0,1)];
    end
    if numGroups >= 12
        customColorMap = [customColorMap; BF_getcmap('purplered',numColorMapGrads,0,1)];
    end
    if numGroups >= 13
        customColorMap = [customColorMap; BF_getcmap('redpurple',numColorMapGrads,0,1)];
    end
    if numGroups >= 14
        customColorMap = [customColorMap; BF_getcmap('orangered',numColorMapGrads,0,1)];
    end
    if numGroups >= 15
        customColorMap = [customColorMap; BF_getcmap('yelloworangebrown',numColorMapGrads,0,1)];
    end
    if numGroups >= 16
        customColorMap = [customColorMap; BF_getcmap('greenblue',numColorMapGrads,0,1)];
    end
    if numGroups >= 17
        customColorMap = [customColorMap; BF_getcmap('bluepurple',numColorMapGrads,0,1)];
    end
    if numGroups >= 18
        customColorMap = [customColorMap; BF_getcmap('bluegreen',numColorMapGrads,0,1)];
    end
    if numGroups >= 19
        fprintf(1,'Too many data groups to colour them correctly\n');
        customColorMap = BF_getcmap('spectral',numColorMapGrads);
    end
end

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
f = figure('color','w');

if addTimeSeries
    %  First make an additional plot to include time-series subsegments

    sp1 = subplot(1,5,1); ax1 = gca; box('on'); hold on
    ax1.YTick = (1:numTS);
    ax1.YTickLabel = {TimeSeries.Name};
    ax1.YLim = [0.5,numTS+0.5];
    ax1.XLim = [1,timeSeriesLength];
    xlabel('Time (samples)');
    ax1.TickLabelInterpreter = 'none';
    NormMinMax = @(x) (x-min(x))/(max(x)-min(x));
    for j = 1:numTS
        tsData = TimeSeries(j).Data(1:timeSeriesLength);
        lengthHere = min(timeSeriesLength,length(tsData));
        plot(1:lengthHere,j-0.5+NormMinMax(tsData),'-k');
        if j < numTS
            plot([1,timeSeriesLength],(j+0.5)*ones(2,1),':k')
        end
    end

    sp2 = subplot(1,5,2:5); ax2 = gca; box('on'); hold on
    linkaxes([ax1,ax2],'y');
end

% ------------------------------------------------------------------------------
%% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
% (alternative is to use imagesc)
colormap(customColorMap)
imagesc(TS_DataMat);
% pcolor([TS_DataMat, zeros(size(TS_DataMat,1),1); zeros(1,size(TS_DataMat,2)+1)]);
% shading flat

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

title(sprintf('Data matrix (%u x %u)',numTS,numOps))

if addTimeSeries
    ax2.YTickLabel = {};

    % Reposition tight
    sp2.Position = [0.4,0.1,0.55,0.85];
    sp1.Position = [0.05,0.1,0.2,0.85];
    % Tight on left:
    sp1.Position(1) = sp2.Position(1) - sp1.Position(3) - 0.01;
    % Both have the same y and height:
    sp1.Position(2) = sp2.Position(2); % y
    sp1.Position(4) = sp2.Position(4); % height
else
    % Rows -- time series need labels
    ylabel('Time series')
    ax2.YTickLabel = {TimeSeries.Name};
end

end
