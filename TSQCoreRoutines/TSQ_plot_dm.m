function TSQ_plot_dm(norcl,kwgs,gi,F,customorder,customcmap)
% Plot the data matrix (time series versus operations) table
% Ben Fulcher 9/12/2009 -- updated for mySQL system
% Ben Fulcher 31/3/2010 -- added norcl option
% Ben Fulcher 24/6/2010 -- cleaned up the axis labels fiasco, fixed error
%                           in which customordering didn't change this when
%                           it should. Redid the colormap handling of
%                           multiple groups.
% Ben Fulcher 25/6/2010 -- added F input

%% OUTPUT
% Gives a colormap plot of time series (rows) and operations (columns)

%% Fill in default values if labelx unspecified; check for valid inputs
% Display as the normalized or clustered data
% What sorts of axes labels are appropriate
if nargin < 1 || isempty(norcl)
    norcl = 'cl'; % by default visualize the clustered matrix
end
if nargin < 2; kwgs = {}; end % no groups

% Differential colouring of keyword groups?
if ischar(kwgs)
    kwgs = {kwgs};
end
if nargin < 3
    gi = []; % automatically get indicies if necessary
end 
if nargin < 4
   F = []; % load from TS_loc_N or TS_loc_cl
end
if nargin < 5 || isempty(customorder)
	customorder = {[],[]};
end
if nargin < 6
    customcmap = 'redyellowblue';
end

%% Read in the data
if isstruct(norcl)
    % can specify all of this in the norcl argument
    tsf = norcl.tsf;
    mlab = norcl.mlab;
    F = norcl.F;
else
    fprintf(1,'Reading data and guides from file...\n');
    if strcmp(norcl,'cl')
        if isempty(F)
            a = which('TS_loc_cl.mat'); % first check it exists
            if isempty(a), error('TS_loc_cl not found: run TSQ_cluster'); end
            load TS_loc_cl.mat TS_loc_cl
            F = TS_loc_cl; clear TS_loc_cl
        end
        load TS_loc_guides_cl.mat tsfcl mlabcl
        tsf = tsfcl; clear tsfcl
        mlab = mlabcl; clear mlabcl
    elseif strcmp(norcl,'norm')
        if isempty(F)
            a = which('TS_loc_N.mat');
            if isempty(a), error('TS_loc_N not found: run TSQ_normalize'); end
            load TS_loc_N.mat TS_loc_N
            F = TS_loc_N; clear TS_loc_N
        end
        load TS_loc_guides_N.mat tsfn mlabn
        tsf = tsfn; clear tsfn
        mlab = mlabn; clear mlabn
    else
        error('Unknown specifier %s, please specify ''norm'' or ''cl''',norcl)
    end
end

[nts, nops] = size(F); % label in this way -- ts as rows

%% Reorder according to customorder
if ~isempty(customorder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified\n');
	F = F(customorder{1},:);
    tsf = tsf(customorder{1});
end

if ~isempty(customorder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified\n');
	F = F(:,customorder{2});
    mlab = mlab(customorder{2});
end


%% Plot the object in a new figure
figure('color','w'); box('on');
title(sprintf('Data matrix of size %u x %u',nts,nops))
ng = 6; % number of gradations in each set of colourmap
nkwgs = size(kwgs,1); % number of keyword groups

if nkwgs > 0
    fprintf(1,'Dividing data into %u groups\n',nkwgs);

    if isempty(gi)
        gi = SUB_autolabelQ(kwgs,metorts,norcl);
    end
    Ng = length(gi);
    
    % add a group for unlabelled data items if they exist
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
    
    %% Change range of F to make use of new colormap appropriately
    ff = 0.9999999;
    squashme = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    F = squashme(F);
    for jo = 1:Ng
        F(gi{jo},:) = squashme(F(gi{jo},:)) + jo - 1;
    end
end
Ng = length(gi);
% set the colormap
if Ng <= 1
    if strcmp(customcmap,'redyellowblue');
        customcmap = BF_getcmap('redyellowblue',ng,0);
    else
        customcmap = gray(ng);
    end
    colormap(customcmap)
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

% have to surround by zeros for an accurate pcolor:
pcolor([F, zeros(size(F,1),1); zeros(1,size(F,2)+1)]);
shading flat

%%% Format the plot
% Axis labels:
set(gca,'YTick',1 + (0.5:1:size(F,1)),'YTickLabel',tsf); % time series
if nops < 1000 % otherwise don't bother
    set(gca,'XTick',1 + (0.5:1:size(F,2)),'XTickLabel',mlab);
end
title(sprintf('Data matrix (%ux%u)',size(F,1),size(F,2)))
set(gca,'FontSize',8) % set font size

end