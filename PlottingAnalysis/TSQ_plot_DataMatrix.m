% ------------------------------------------------------------------------------
% TSQ_plot_DataMatrix
% ------------------------------------------------------------------------------
% 
% Plot the data matrix.
% 
%---INPUTS:
% WhatData: specify 'norm' for normalized data in HCTSA_N.mat, 'cl' for clustered
%         data in HCTSA_cl.mat (default)
% kwgs:  specify keyword groups to color different keywords differently in the
%         data matrix [opt].
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

function TSQ_plot_DataMatrix(WhatData,ColorGroups,CustomOrder,CustomColorMap)

% ------------------------------------------------------------------------------
%% Check Inputs:
% ------------------------------------------------------------------------------
% What data to plot:
if nargin < 1 || isempty(WhatData)
    WhatData = 'cl'; % Load data from HCTSA_cl by default
end
% Color groups of time series differently:
if nargin < 2 || isempty(ColorGroups)
    ColorGroups = 0; % Don't color groups by default
end
if nargin < 3 || isempty(CustomOrder)
	CustomOrder = {[],[]};
end
if nargin < 4
    CustomColorMap = 'redyellowblue';
end

% --------------------------------------------------------------------------
%% Read in the data
% --------------------------------------------------------------------------
if isstruct(WhatData)
    % can specify all of this in the WhatData argument
    TimeSeries = WhatData.TimeSeries;
    Operations = WhatData.Operations;
    TS_DataMat = WhatData.TS_DataMat;
else
    if strcmp(WhatData,'cl')
        TheFile = 'HCTSA_cl.mat'; TheRoutine = 'TSQ_cluster';
    elseif strcmp(WhatData,'norm')
        TheFile = 'HCTSA_N.mat'; TheRoutine = 'TSQ_normalize';
    elseif ischar(WhatData) % Specify a filename
        a = which(TheFile); % First check it exists
        if isempty(a)
            error('\n%s not found. You should probably run %s...',TheFile,TheRoutine);
        end
    else
        error(['Unknown specifier ''%s'', please input the data structure, ' ...
                        '\t\nor specify ''norm'', ''cl'', or a FileName'],WhatData)
    end
    fprintf(1,'Reading data from %s...',TheFile);
    load(TheFile,'TimeSeries','Operations','TS_DataMat')
    fprintf(1,' Done.\n');
end

if ColorGroups
    TimeSeriesGroups = [TimeSeries.Group];
    fprintf(1,'Coloring groups of time series...\n');
end

TimeSeriesFileNames = {TimeSeries.FileName}; clear TimeSeries; % Just extract filenames
OperationNames = {Operations.Name}; clear Operations; % Just extract operation names

[nts, nops] = size(TS_DataMat); % size of the data matrix

% ------------------------------------------------------------------------------
%% Reorder according to CustomOrder
% ------------------------------------------------------------------------------
if ~isempty(CustomOrder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified\n');
	TS_DataMat = TS_DataMat(CustomOrder{1},:);
    TimeSeriesFileNames = TimeSeriesFileNames(CustomOrder{1});
end

if ~isempty(CustomOrder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified\n');
	TS_DataMat = TS_DataMat(:,CustomOrder{2});
    OperationNames = OperationNames(CustomOrder{2});
end

% --------------------------------------------------------------------------
%% Plot the data matrix in a new figure
% --------------------------------------------------------------------------
figure('color','w'); box('on');
title(sprintf('Data matrix of size %u x %u',nts,nops))
ng = 6; % number of gradations in each set of colourmap

if ColorGroups    
    gi = BF_ToGroup(TimeSeriesGroups);
    
    NumGroups = length(gi);
    
    % Add a group for unlabelled data items if they exist
    if sum(cellfun(@length,gi)) < nts
        % Add an unlabelled class
        gi0 = gi;
        gi = cell(NumGroups+1,1);
        for i = 1:NumGroups
            gi{i} = gi0{i};
        end
        clear gi0;
        gi{end} = setxor(1:nts,vertcat(gi{:}));
        NumGroups = NumGroups + 1;
    end
    
    fprintf(1,'Coloring data according to %u groups\n',NumGroups);
    
    % Change range of TS_DataMat to make use of new colormap appropriately
    ff = 0.9999999;
    squashme = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    TS_DataMat = squashme(TS_DataMat);
    for jo = 1:NumGroups
        TS_DataMat(gi{jo},:) = squashme(TS_DataMat(gi{jo},:)) + jo - 1;
    end
else
    NumGroups = 0;
end

% --------------------------------------------------------------------------
%% Set the colormap
% --------------------------------------------------------------------------
if NumGroups <= 1
    if strcmp(CustomColorMap,'redyellowblue');
        CustomColorMap = BF_getcmap('redyellowblue',ng,0);
    else
        CustomColorMap = gray(ng);
    end
    colormap(CustomColorMap)
else
    cmap = colormap(BF_getcmap('blues',ng,0,1));
    if NumGroups >= 2
        cmap = [cmap; BF_getcmap('greens',ng,0,1)];
    end
    if NumGroups >= 3
        cmap = [cmap; BF_getcmap('oranges',ng,0,1)];
    end
    if NumGroups >= 4
        cmap = [cmap; BF_getcmap('purples',ng,0,1)];
    end
    if NumGroups >= 5
        cmap = [cmap; BF_getcmap('reds',ng,0,1)];
    end
    if NumGroups >= 6
        cmap = [cmap; pink(ng)];
    end
    if NumGroups >= 7
        cmap = [cmap; gray(ng)];
    end
    if NumGroups >= 8
        cmap = [cmap; BF_getcmap('yelloworangered',ng,0,1)];
    end
    if NumGroups >= 9
        cmap = [cmap; BF_getcmap('purplebluegreen',ng,0,1)];
    end
    if NumGroups >= 10
        cmap = [cmap; BF_getcmap('yellowgreenblue',ng,0,1)];
    end
    if NumGroups >= 11
        cmap = [cmap; BF_getcmap('purpleblue',ng,0,1)];
    end
    if NumGroups >= 12
        cmap = [cmap; BF_getcmap('purplered',ng,0,1)];
    end
    if NumGroups >= 13
        cmap = [cmap; BF_getcmap('redpurple',ng,0,1)];
    end
    if NumGroups >= 14
        cmap = [cmap; BF_getcmap('orangered',ng,0,1)];
    end
    if NumGroups >= 15
        cmap = [cmap; BF_getcmap('yelloworangebrown',ng,0,1)];
    end
    if NumGroups >= 16
        cmap = [cmap; BF_getcmap('greenblue',ng,0,1)];
    end
    if NumGroups >= 17
        cmap = [cmap; BF_getcmap('bluepurple',ng,0,1)];
    end
    if NumGroups >= 18
        cmap = [cmap; BF_getcmap('bluegreen',ng,0,1)];
    end
    if NumGroups >= 19
        fprintf(1,'Too many data groups to colour them correctly\n');
        cmap = BF_getcmap('spectral',ng);
    end
    colormap(cmap)
end

% ------------------------------------------------------------------------------
%% Plot the data matrix
% ------------------------------------------------------------------------------
% Surround by zeros for an accurate and inclusive pcolor:
pcolor([TS_DataMat, zeros(size(TS_DataMat,1),1); zeros(1,size(TS_DataMat,2)+1)]);
shading flat

% ------------------------------------------------------------------------------
% Superimpose green/yellow rectangles over NaN values
% ------------------------------------------------------------------------------
if any(isnan(TS_DataMat(:)))
    Green = BF_getcmap('greens',5,1);
    [theNaNs_i,theNaNs_j] = find(isnan(TS_DataMat));
    fprintf(1,['Superimposing green/yellow rectangles over all %u NaNs in ' ...
                                'the data matrix\n'],length(theNaNs_i));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i),theNaNs_i(i),1,1],'FaceColor',Green{end}, ...
                        'EdgeColor','y')
    end
end

% --------------------------------------------------------------------------
%% Format the plot
% --------------------------------------------------------------------------
% Axis labels:
set(gca,'YTick',1 + (0.5:1:size(TS_DataMat,1)),'YTickLabel',TimeSeriesFileNames); % time series
if nops < 1000 % otherwise don't bother
    set(gca,'XTick',1 + (0.5:1:size(TS_DataMat,2)),'XTickLabel',OperationNames);
end
title(sprintf('Data matrix (%ux%u)',size(TS_DataMat,1),size(TS_DataMat,2)))
set(gca,'FontSize',8) % set font size

end