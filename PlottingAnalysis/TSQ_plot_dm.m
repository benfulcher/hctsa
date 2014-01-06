% TSQ_plot_dm
% 
% Plot the data matrix.
% 
%---INPUTS:
% norcl: specify 'norm' for normalized data in HCTSA_N.mat, 'cl' for clustered
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
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TSQ_plot_dm(norcl,ColorGroups,TS_DataMat,CustomOrder,CustomColorMap)

% Visualize normalized or clustered matrix (clustered by default)
if nargin < 1 || isempty(norcl)
    norcl = 'cl';
end
if nargin < 2 || isempty(ColorGroups)
    ColorGroups = 0; % don't color groups
end
% Differential colouring of keyword groups?
% if nargin < 2; kwgs = {}; end % no groups
% if ischar(kwgs)
%     kwgs = {kwgs};
% end
% if nargin < 3
%     gi = []; % automatically get indicies if necessary
% end 
if nargin < 4
   TS_DataMat = []; % load from TS_loc_N or TS_loc_cl
end
if nargin < 5 || isempty(CustomOrder)
	CustomOrder = {[],[]};
end
if nargin < 6
    CustomColorMap = 'redyellowblue';
end

%% Read in the data
if isstruct(norcl)
    % can specify all of this in the norcl argument
    TimeSeriesFileNames = norcl.TimeSeriesFileNames;
    OperationNames = norcl.OperationNames;
    TS_DataMat = norcl.F;
else
    if strcmp(norcl,'cl')
        TheFile = 'HCTSA_cl.mat'; TheRoutine = 'TSQ_cluster';
    elseif strcmp(norcl,'norm')
        TheFile = 'HCTSA_N.mat'; TheRoutine = 'TSQ_normalize';
    else
        error('Unknown specifier %s, please specify ''norm'' or ''cl''',norcl)
    end
    
    fprintf(1,'Reading data from %s...',TheFile);
    
    if isempty(TS_DataMat)
        a = which(TheFile); % First check it exists
        if isempty(a)
            error('\n%s not found. You should probably run %s...',TheFile,TheRoutine);
        end
        load(TheFile,'TS_DataMat')
        fprintf(1,' Done.\n');
    end
    load(TheFile,'TimeSeries','Operations')
end
if ColorGroups
    load(TheFile,'GroupNames')
    TimeSeriesGroups = [TimeSeries.Group];
end
TimeSeriesFileNames = {TimeSeries.FileName}; clear TimeSeries; % Just extract filenames
OperationNames = {Operations.Name}; clear Operations; % Just extract operation names


[nts, nops] = size(TS_DataMat); % label in this way -- ts as rows

%% Reorder according to CustomOrder
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

%% Plot the object in a new figure
figure('color','w'); box('on');
title(sprintf('Data matrix of size %u x %u',nts,nops))
ng = 6; % number of gradations in each set of colourmap

if ColorGroups
    Ng = length(GroupNames); % number of keyword groups
    
    fprintf(1,'Coloring data according to %u groups\n',Ng);
    
    gi = BF_ToGroup(TimeSeriesGroups);
    
    % Add a group for unlabelled data items if they exist
    if sum(cellfun(@length,gi)) < nts
        % we need to add an unlabelled class
        gi0 = gi;
        gi = cell(Ng+1,1);
        for i = 1:Ng
            gi{i} = gi0{i};
        end
        clear gi0;
        gi{end} = setxor(1:nts,cell2mat(gi));
        Ng = Ng + 1;
    end
    
    %% Change range of TS_DataMat to make use of new colormap appropriately
    ff = 0.9999999;
    squashme = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    TS_DataMat = squashme(TS_DataMat);
    for jo = 1:Ng
        TS_DataMat(gi{jo},:) = squashme(TS_DataMat(gi{jo},:)) + jo - 1;
    end
end
Ng = length(gi);
% set the colormap
if Ng <= 1
    if strcmp(CustomColorMap,'redyellowblue');
        CustomColorMap = BF_getcmap('redyellowblue',ng,0);
    else
        CustomColorMap = gray(ng);
    end
    colormap(CustomColorMap)
else
    cmap = colormap(BF_getcmap('blues',ng,0,1));
    if Ng >= 2
        cmap = [cmap; BF_getcmap('greens',ng,0,1)];
    end
    if Ng >= 3
        cmap = [cmap; BF_getcmap('oranges',ng,0,1)];
    end
    if Ng >= 4
        cmap = [cmap; BF_getcmap('purples',ng,0,1)];
    end
    if Ng >= 5
        cmap = [cmap; BF_getcmap('reds',ng,0,1)];
    end
    if Ng >= 6
        cmap = [cmap;pink(ng)];
    end
    if Ng >= 7
        cmap = [cmap;gray(ng)];
    end
    if Ng >= 8
        cmap = [cmap; BF_getcmap('yelloworangered',ng,0,1)];
    end
    if Ng >= 9
        cmap = [cmap; BF_getcmap('purplebluegreen',ng,0,1)];
    end
    if Ng >= 10
        cmap = [cmap; BF_getcmap('yellowgreenblue',ng,0,1)];
    end
    if Ng >= 11
        cmap = [cmap; BF_getcmap('purpleblue',ng,0,1)];
    end
    if Ng >= 12
        cmap = [cmap; BF_getcmap('purplered',ng,0,1)];
    end
    if Ng >= 13
        cmap = [cmap; BF_getcmap('redpurple',ng,0,1)];
    end
    if Ng >= 14
        cmap = [cmap; BF_getcmap('orangered',ng,0,1)];
    end
    if Ng >= 15
        cmap = [cmap; BF_getcmap('yelloworangebrown',ng,0,1)];
    end
    if Ng >= 16
        cmap = [cmap; BF_getcmap('greenblue',ng,0,1)];
    end
    if Ng >= 17
        cmap = [cmap; BF_getcmap('bluepurple',ng,0,1)];
    end
    if Ng >= 18
        cmap = [cmap; BF_getcmap('bluegreen',ng,0,1)];
    end
    if Ng >= 19
        fprintf(1,'Too many data groups to colour them correctly\n');
        cmap = BF_getcmap('spectral',ng);
    end
    colormap(cmap)
end

% Surround by zeros for an accurate pcolor:
pcolor([TS_DataMat, zeros(size(TS_DataMat,1),1); zeros(1,size(TS_DataMat,2)+1)]);
shading flat

%%% Format the plot
% Axis labels:
set(gca,'YTick',1 + (0.5:1:size(TS_DataMat,1)),'YTickLabel',TimeSeriesFileNames); % time series
if nops < 1000 % otherwise don't bother
    set(gca,'XTick',1 + (0.5:1:size(TS_DataMat,2)),'XTickLabel',OperationNames);
end
title(sprintf('Data matrix (%ux%u)',size(TS_DataMat,1),size(TS_DataMat,2)))
set(gca,'FontSize',8) % set font size

end