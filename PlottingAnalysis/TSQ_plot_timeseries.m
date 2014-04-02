% ------------------------------------------------------------------------------
% TSQ_plot_timeseries
% ------------------------------------------------------------------------------
% 
% Plots the time series read from a local file, in a specified format.
% 
%---INPUTS:
% WhatData, The data to get information from: can be a structure, or 'norm' or
%           'cl' to load from HCTSA_N or HCTSA_cl
% WhatTimeSeries, Can provide indices to plot that subset, a keyword to plot
%                   matches to the keyword, 'all' to plot all, or an empty vector
%                   to plot default groups in TimeSeries.Group
% NumPerGroup, If plotting groups, plots this many examples per group
% maxL, the maximum number of samples of each time series to plot
% titleme, shows time-series labels on the plot.
% opts, additional plotting options as a structure.
% 
%----HISTORY:
% Previously called 'TSQ_plot_examples'
% Ben Fulcher, 9/4/2010
% Ben Fulcher, 13/5/2010 added F option (a matrix, 'norm', or 'cl')
% Ben Fulcher, 24/6/2010 added option to show examples from each class in
%                       kwgs, rather than all from the first class. In this
%                       case, <NumPerGroup> means per class.
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

function TSQ_plot_timeseries(WhatData,WhatTimeSeries,NumPerGroup,maxL,titleme,opts)

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
% F -- get guide from 'norm' or 'cl'
if nargin < 1 || isempty(WhatData)
    WhatData = 'norm';
end

% Can specify a reduced set of time series by keyword
if nargin < 2
    WhatTimeSeries = '';
end

if nargin < 3 || isempty(NumPerGroup)
    % Default: plot 5 time series per group
    NumPerGroup = 5;
end

if nargin < 4
    % Maximum length of time series to display (otherwise crops)
    % If empty, displays all of all time series
    maxL = [];
end

if nargin < 5 || isempty(titleme)
    % Show time series labels -- allows more to be fit into plot
    titleme = 1;
end

if nargin < 6
	opts = [];
end

% ------------------------------------------------------------------------------
% Evaluate any custom plotting options specified in the structure opts
% ------------------------------------------------------------------------------

if isstruct(opts) && isfield(opts,'howtofilter')
    howtofilter = opts.howtofilter;
else
    howtofilter = 'rand'; % by default
end
if isstruct(opts) && isfield(opts,'gic')
    gic = opts.gic; % local color labels -- vector
else
    gic = [];
end
% Specify the colormap to use
if isstruct(opts) && isfield(opts,'cmap')
    cmap = opts.cmap;
else
    cmap = 'set1';
end
% Specify whether to make a free-form plot
if isstruct(opts) && isfield(opts,'freeform')
    freeform = opts.freeform;
else
    freeform = 0; % do a normal subplotted figure
end
% Specify line width for plotting
if isstruct(opts) && isfield(opts,'LineWidth')
    lw = opts.LineWidth;
else
    lw = 1; % do a normal subplotted figure
end

% ------------------------------------------------------------------------------
%% Load data
% ------------------------------------------------------------------------------
if isstruct(WhatData)
    % Provide it all yourself
    TimeSeries = WhatData.TimeSeries;
else    
    if strcmp(WhatData,'cl')
        TheFile = 'HCTSA_cl.mat';
    elseif strcmp(WhatData,'norm')
        TheFile = 'HCTSA_N.mat';
    end
    load(TheFile,'TimeSeries');
end


% ------------------------------------------------------------------------------
%% Get group indices:
% ------------------------------------------------------------------------------
if isempty(WhatTimeSeries)
    % Use default groups
    GroupIndices = BF_ToGroup([TimeSeries.Group]);
    fprintf(1,'Plotting from %u groups of time series from file\n',length(GroupIndices));
elseif strcmp(WhatTimeSeries,'all')
    % Plot all time series
    GroupIndices = {1:length(TimeSeries)};
elseif ischar(WhatTimeSeries)
    % Just plot this group
    % First load group names:
    if isstruct(WhatData)
        GroupNames = WhatData.GroupNames;
    else
        load(TheFile,'GroupNames');
    end
    a = strcmp(WhatTimeSeries,GroupNames);
    GroupIndices = {find([TimeSeries.Group]==find(a))};
    fprintf(1,'Plotting %u time series matching group name ''%s''\n',length(GroupIndices{1}),WhatTimeSeries);
else % Provided a custom range as a vector
    GroupIndices = {WhatTimeSeries};
    fprintf(1,'Plotting the %u time series matching indices provided\n',length(WhatTimeSeries));
end
NumGroups = length(GroupIndices);

% ------------------------------------------------------------------------------
%% Do the plotting
% ------------------------------------------------------------------------------
% Want NumPerGroup from each time series group
iplot = zeros(NumGroups*NumPerGroup,1);
classes = zeros(NumGroups*NumPerGroup,1);
nhere = zeros(NumGroups,1);
gil = cellfun(@length,GroupIndices);
% howtofilter = 'rand';
for i = 1:NumGroups
    % filter down to NumPerGroup if too many in group, otherwise plot all in
    % group
    switch howtofilter
        
        case 'firstcome'
            % just plot first in group (useful when ordered by closeness to
            % cluster centre)
            jj = (1:min(NumPerGroup,gil(i)));
            
        case 'rand'
            % select ones to plot at random
            if gil(i) > NumPerGroup
                jj = randperm(gil(i)); % randomly selected
                if length(jj) > NumPerGroup
                    jj = jj(1:NumPerGroup);
                end
            else
                jj = (1:min(NumPerGroup,gil(i))); % retain order if not subsampling
            end
    end
    nhere(i) = length(jj); % could be less than NumPerGroup if a smaller group
    rh = sum(nhere(1:i-1))+1:sum(nhere(1:i)); % range here
    iplot(rh) = GroupIndices{i}(jj);
    classes(rh) = i;
end
rkeep = (iplot > 0);
classes = classes(rkeep);
iplot = iplot(rkeep); % contains all the indicies of time series to plot (in order)
Nplot = length(iplot);


fprintf(1,'Plotting %u (/%u) time series from %u classes\n', ...
                    Nplot,sum(cellfun(@length,GroupIndices)),NumGroups)
if ~isempty(gic)
    classes = gic; % override with our group information
    % better be firstcome, and makes sense that GroupIndices is just a vector of
    % indicies, otherwise group information is already available!
    if iscell(cmap)
        c = cmap; % specified custom colors as a cell
    else
        c = BF_getcmap(cmap,max(classes),1); % 'set2'
    end
else
    if NumGroups==1;
        c = {'k'}; % just plot in black
    else
        if iscell(cmap)
            c = cmap;
        else
            c = BF_getcmap(cmap,NumGroups,1); % 'set2'
        end
    end    
end

figure('color','w'); box('on');
Ls = zeros(Nplot,1);
if freeform
	% FREEFORM: make all within a single plot with text labels
    hold on;
	yr = linspace(1,0,Nplot+1);
    inc = abs(yr(2)-yr(1)); % size of increment
    yr = yr(2:end);
	ls = zeros(Nplot,1); % lengths of each time series
	if isempty(maxL)
		for i = 1:Nplot
			ls(i) = length(TimeSeries(iplot(i)).Data);
		end
		maxN = max(ls); % maximum length of all time series to plot
	else
		maxN = maxL;
	end
	% Set up axes ticks
	set(gca,'XTick',linspace(0,1,3),'XTickLabel',round(linspace(0,maxN,3)))
	set(gca,'YTick',[],'YTickLabel',{})
	
	for i = 1:Nplot
	    fn = TimeSeries(iplot(i)).FileName; % the filename
	    kw = TimeSeries(iplot(i)).Keywords; % the keywords
	    x = TimeSeries(iplot(i)).Data;
	    N0 = length(x);
		% rectangle('Position',[0,yr(i),1,inc],'EdgeColor','k')
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
		
		plot(xx,xsc,'-','color',c{classes(i)},'LineWidth',lw)

	    % annotate text labels
		if titleme
			thetit = sprintf('{%u} %s [%s] (%u)',TimeSeries(iplot(i)).ID,fn,kw,N0);
			text(0.01,yr(i)+0.9*inc,thetit,'interpreter','none','FontSize',8)
	    end
	end
	xlim([0,1]) % Don't let the axes annoyingly slip out
	
else
	for i = 1:Nplot
	    subplot(Nplot,1,i)
	    fn = TimeSeries(iplot(i)).FileName; % the filename
	    kw = TimeSeries(iplot(i)).Keywords; % the keywords
	    x = TimeSeries(iplot(i)).Data;
	    N = length(x);
    
	    % prepare text for the title
		if titleme
			startbit = sprintf('{%u} %s [%s]',TimeSeries(iplot(i)).ID,fn,kw);
	    end

	    % Plot the time series
	    if isempty(maxL)
	        % no maximum length specified: plot the whole time series
	        plot(x,'-','color',c{classes(i)})
	        Ls(i) = N;
	        if titleme
	            title([startbit ' (' num2str(N) ')'],'interpreter','none','FontSize',8);
	        end
	%         set(gca,'Xtick',[],'FontSize',8) % don't put ticks on any time series
	    else
	        % Specified a maximum length of time series to plot: maxL
	        if N <= maxL
	            plot(x,'-','color',c{classes(i)});
	            Ls(i) = N;
	            if titleme
	                title([startbit ' (' num2str(N) ')'],'interpreter','none','FontSize',8);
	            end
	        else
	            sti = randi(N-maxL,1);
	            plot(x(sti:sti+maxL),'-','color',c{classes(i)}) % plot a random maxL-length portion of the time series
	            Ls(i) = maxL;
	            if titleme
	                title([startbit ' (' num2str(N) ' :: ' num2str(sti) '-' num2str(sti+maxL) ')'],'interpreter','none','FontSize',8);
	            end
	        end
	    end
	    set(gca,'YTickLabel','');
	    if i~=Nplot
	        set(gca,'XTickLabel','','FontSize',8) % put the ticks for the last time series
	    end
	%     set(gca,'Ytick',[])
	end
	
	
	% Set all xlims so that they have the same x-axis limits
	for i = 1:Nplot
	    subplot(Nplot,1,i); set(gca,'xlim',[1,max(Ls)])
	end
end


end