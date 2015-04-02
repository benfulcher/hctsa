% ------------------------------------------------------------------------------
% TSQ_plot_DataMatrix
% ------------------------------------------------------------------------------
% 
% Plot the data matrix.
% 
%---INPUTS:
% whatData: specify 'norm' for normalized data in HCTSA_N.mat, 'cl' for clustered
%         data in HCTSA_cl.mat (default), or specify a filename to load data
%         from that file.
% colorGroups: Set to 1 to color time-series groups with different colormaps.
% customOrder: reorder rows and columns according to provided permutation vectors
% customColorMap: use a custom color map (name to match an option in BF_getcmap)
% colorNaNs: whether to plot rectangles over special-values in the matrix (default: 1)
% 
%---OUTPUT:
% Produces a colormap plot of the data matrix with time series as rows and
%   operations as columns.
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

function TSQ_plot_DataMatrix(whatData,colorGroups,customOrder,customColorMap,colorNaNs)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
% What data to plot:
if nargin < 1 || isempty(whatData)
    whatData = 'cl'; % Load data from HCTSA_cl by default
end
% Color groups of time series differently:
if nargin < 2 || isempty(colorGroups)
    colorGroups = 0; % Don't color groups by default
end
if nargin < 3 || isempty(customOrder)
	customOrder = {[],[]};
end
if nargin < 4 || isempty(customColorMap)
    customColorMap = 'redyellowblue';
end
if nargin < 5 || isempty(colorNaNs)
    colorNaNs = 1;
end

% --------------------------------------------------------------------------
%% Read in the data
% --------------------------------------------------------------------------
if isstruct(whatData)
    % Can specify all of these fields in the whatData argument
    TimeSeries = whatData.TimeSeries;
    Operations = whatData.Operations;
    TS_DataMat = whatData.TS_DataMat;
else
    if strcmp(whatData,'cl')
        theFile = 'HCTSA_cl.mat'; TheRoutine = 'TSQ_cluster';
    elseif strcmp(whatData,'norm')
        theFile = 'HCTSA_N.mat'; TheRoutine = 'TSQ_normalize';
    elseif ischar(whatData) && exist(whatData) % Specify a filename
        theFile = whatData;
        a = which(theFile); % First check it exists
        if isempty(a)
            error('\n%s not found. You should probably run %s...',theFile,TheRoutine);
        end
    else
        error(['Unknown specifier ''%s'', please input the data structure, ' ...
                        '\t\n or specify ''norm'', ''cl'', or a FileName'],whatData)
    end
    fprintf(1,'Reading data from %s...',theFile);
    load(theFile,'TimeSeries','Operations','TS_DataMat')
    fprintf(1,' Done.\n');
end

if colorGroups
    timeSeriesGroups = [TimeSeries.Group];
    fprintf(1,'Coloring groups of time series...\n');
end

timeSeriesFileNames = {TimeSeries.FileName}; clear TimeSeries; % Just extract filenames
operationNames = {Operations.Name}; clear Operations; % Just extract operation names

[numTS, numOps] = size(TS_DataMat); % size of the data matrix

% ------------------------------------------------------------------------------
%% Reorder according to customOrder
% ------------------------------------------------------------------------------
if ~isempty(customOrder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified.\n');
	TS_DataMat = TS_DataMat(customOrder{1},:);
    timeSeriesFileNames = timeSeriesFileNames(customOrder{1});
end

if ~isempty(customOrder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified.\n');
	TS_DataMat = TS_DataMat(:,customOrder{2});
    operationNames = operationNames(customOrder{2});
end

% --------------------------------------------------------------------------
%% Plot the data matrix in a new figure
% --------------------------------------------------------------------------
figure('color','w'); box('on');
title(sprintf('Data matrix of size %u x %u',numTS,numOps))
ng = 6; % number of gradations in each set of colourmap

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
    squashme = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    TS_DataMat = squashme(TS_DataMat);
    for jo = 1:numGroups
        TS_DataMat(gi{jo},:) = squashme(TS_DataMat(gi{jo},:)) + jo - 1;
    end
else
    numGroups = 0;
end

% --------------------------------------------------------------------------
%% Set the colormap
% --------------------------------------------------------------------------
if numGroups <= 1
    if strcmp(customColorMap,'redyellowblue');
        customColorMap = BF_getcmap('redyellowblue',ng,0);
    else
        customColorMap = gray(ng);
    end
    colormap(customColorMap)
else
    cmap = colormap(BF_getcmap('blues',ng,0,1));
    if numGroups >= 2
        cmap = [cmap; BF_getcmap('greens',ng,0,1)];
    end
    if numGroups >= 3
        cmap = [cmap; BF_getcmap('oranges',ng,0,1)];
    end
    if numGroups >= 4
        cmap = [cmap; BF_getcmap('purples',ng,0,1)];
    end
    if numGroups >= 5
        cmap = [cmap; BF_getcmap('reds',ng,0,1)];
    end
    if numGroups >= 6
        cmap = [cmap; pink(ng)];
    end
    if numGroups >= 7
        cmap = [cmap; gray(ng)];
    end
    if numGroups >= 8
        cmap = [cmap; BF_getcmap('yelloworangered',ng,0,1)];
    end
    if numGroups >= 9
        cmap = [cmap; BF_getcmap('purplebluegreen',ng,0,1)];
    end
    if numGroups >= 10
        cmap = [cmap; BF_getcmap('yellowgreenblue',ng,0,1)];
    end
    if numGroups >= 11
        cmap = [cmap; BF_getcmap('purpleblue',ng,0,1)];
    end
    if numGroups >= 12
        cmap = [cmap; BF_getcmap('purplered',ng,0,1)];
    end
    if numGroups >= 13
        cmap = [cmap; BF_getcmap('redpurple',ng,0,1)];
    end
    if numGroups >= 14
        cmap = [cmap; BF_getcmap('orangered',ng,0,1)];
    end
    if numGroups >= 15
        cmap = [cmap; BF_getcmap('yelloworangebrown',ng,0,1)];
    end
    if numGroups >= 16
        cmap = [cmap; BF_getcmap('greenblue',ng,0,1)];
    end
    if numGroups >= 17
        cmap = [cmap; BF_getcmap('bluepurple',ng,0,1)];
    end
    if numGroups >= 18
        cmap = [cmap; BF_getcmap('bluegreen',ng,0,1)];
    end
    if numGroups >= 19
        fprintf(1,'Too many data groups to colour them correctly\n');
        cmap = BF_getcmap('spectral',ng);
    end
    colormap(cmap)
end

% ------------------------------------------------------------------------------
%% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
% (alternative is to use imagesc)
pcolor([TS_DataMat, zeros(size(TS_DataMat,1),1); zeros(1,size(TS_DataMat,2)+1)]);
shading flat

% ------------------------------------------------------------------------------
% Superimpose colored rectangles over NaN values
% ------------------------------------------------------------------------------
if colorNaNs && any(isnan(TS_DataMat(:)))
    [theNaNs_i,theNaNs_j] = find(isnan(TS_DataMat));
    fprintf(1,['Superimposing black rectangles over all %u NaNs in ' ...
                                'the data matrix\n'],length(theNaNs_i));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i),theNaNs_i(i),1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end
end

% --------------------------------------------------------------------------
%% Format the plot
% --------------------------------------------------------------------------
% Axis labels:
set(gca,'YTick',1 + (0.5:1:size(TS_DataMat,1)),'YTickLabel',timeSeriesFileNames); % time series
if numOps < 1000 % otherwise don't bother
    set(gca,'XTick',1 + (0.5:1:size(TS_DataMat,2)),'XTickLabel',operationNames);
end
title(sprintf('Data matrix (%ux%u)',size(TS_DataMat,1),size(TS_DataMat,2)))
set(gca,'FontSize',8) % set font size
set(gca,'TickLabelInterpreter','none') % Stop from displaying underscores as subscripts
xlabel('Operations')
ylabel('Time series')

end