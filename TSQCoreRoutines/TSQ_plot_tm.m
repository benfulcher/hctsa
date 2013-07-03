function TSQ_plot_tm(kwgs,gi,metorts,labelm,labelx,customorder,norcl,F,customcmap)
%%% TS_plot_tm
% Plot the (time series versus metric) table
%% INPUTS
% 1) kwgs -- keyword groups in which to differentially colour the
%               correlation matrix
% 2) gi -- the indicies of kwgs in the clustered label system (from SUB_autolabelQ, e.g.)
% 3) metorts -- whether the keywords are of time series or operations
% 4) labelm -- method with which to show labels on the plot: 
% 				* 'simple' (default) -- just shows the filenames/metric names
% 				* 'ckey' -- clusters the keywords
% 				* 'skey' -- shows only specific key words
% 5) [opt] labelx -- extra inputs for the label method -- form depends on labelm
% 	'simple' -- should give the (integral) frequency of displaying labels (default is to show all; labelx=1)
% 	'ckey' -- no additional input required
% 	'skey' -- labelx should be a 2-component cell
% 				labelx{1} should be a cell of
% 6) customorder -- specify a custom ordering as {[<row_order>],[<column_order>]}; both given
% 					in the clustered index system; output of TS_prune is ideal for this sort of thing...
% 7) norcl -- specify either 'cl' to visualize the clustered, or 'norm' for
%              the normalized matrix. This could be useful to manually
%              implement a customorder, e.g., or just to see the
%              unclustered data easily

% Ben Fulcher 9/12/2009 -- updated for new mySQL system
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
% Colouring of groups?
% What sorts of axes labels are appropriate
if nargin < 1; kwgs = {}; end % no groups
if nargin < 2; gi = []; end % automatically get indicies if necessary
if nargin < 3 || isempty(metorts); metorts = 'ts'; end % default to time series
if nargin < 4 || isempty(labelm), labelm = 'simple'; end % simple labeling
if nargin < 5 || isempty(labelx) % set labelx
	switch labelm
		case 'simple'
			labelx = 1; % show every label
		case 'ckey'
			% doesn't require any additional inputs
		case 'skey'
			if strcmp(metorts,'timseries'), labelx={{'periodic','noise'},0};
			else labelx = {{'entropy'},0};
			end
			disp('You should really specify inputs for ''skey''');
		otherwise
			disp('No labels'), labelm = 'none';
	end
end
if nargin < 6 || isempty(customorder)
	customorder={[],[]};
end
if nargin < 7 || isempty(norcl)
    norcl = 'cl'; % by default visualize the clustered matrix
end
if nargin < 8
   F = []; % load from TS_loc_N or TS_loc_cl
end
if nargin < 9,
    customcmap = 'redyellowblue';
end
if ~ismember(metorts,{'ops','ts'});
    disp('Choose either ''ops'' or ''ts''');
end
if ischar(kwgs)
    kwgs={kwgs};
end

%% Read in the clustered data
fprintf(1,'Reading data and guides from file...\n');

if isstruct(norcl)
    % can specify all of this in the norcl argument
    tsf = norcl.tsf;
    mlab = norcl.mlab;
    F = norcl.F;
elseif strcmp(norcl,'cl')
    if isempty(F)
        load TS_loc_cl.mat TS_loc_cl
        F = TS_loc_cl; clear TS_loc_cl
    end
    load TS_loc_guides_cl.mat tsfcl mlabcl
    tsf = tsfcl; clear tsfcl
    mlab = mlabcl; clear mlabcl
elseif strcmp(norcl,'norm')
    if isempty(F)
        load TS_loc_N.mat TS_loc_N
        F = TS_loc_N; clear TS_loc_N
    end
    load TS_loc_guides_N.mat tsfn mlabn
    tsf = tsfn; clear tsfn
    mlab = mlabn; clear mlabn
end


if strcmp(metorts,'ops')
    F = F'; % make operations the rows
end

[nts, nm] = size(F); % label in this way -- by default think of ts as rows

%% Reorder according to customorder
if ~isempty(customorder{1}) % reorder rows
	disp('Reordering rows according to custom order');
	F = F(customorder{1},:);
    tsf = tsf(customorder{1});
end

if ~isempty(customorder{2}) % reorder columns
	disp('Reordering columns according to custom order');
	F = F(:,customorder{2});
    mlab = mlab(customorder{2});
end


%% Plot the object in a new figure
figure('color','w'); box('on');
title(sprintf('Data matrix of size %u x %u',nts,nm))
ng = 6; % number of gradations in each set of colourmap
nkwgs = size(kwgs,1); % number of keyword groups

if nkwgs>0
    disp('Looking for; filling keyword groups')
    disp(['Dividing into ' num2str(nkwgs) ' groups']);

    if isempty(gi)
        gi = SUB_autolabelQ(kwgs,metorts,norcl);
    end
    Ng = length(gi);
    
    % if a unlabelled proportion -- add this
    if sum(cellfun(@length,gi)) < nts
        % we need to add an unlabelled class
        gi0 = gi;
        gi = cell(Ng+1,1);
        for i = 1:Ng
            gi{i} = gi0{i};
        end
        clear gi0;
        gi{end} = setxor(1:nts,cell2mat(gi));
        Ng = Ng+1;
    end

    
    %% Change range of F to make use of new colormap appropriately
    ff = 0.9999999;
    squashme = @(x)ff*(x - min(x(:)))/(max(x(:))-min(x(:)));
    F = squashme(F); % (F-min(F(:)))/max(F(:)); % assuming positive values; F now ranges from 0 to ff
%     if strcmp(metorts,'ts');
%     else
%         for jo=1:nkwgs
%             shmoo=F(:,gi{jo});
%             F(:,gi{jo})=(shmoo-min(shmoo(:)))/max(shmoo(:))+jo+1E-6; % also needs to span unit interval
%         end
%     end
    for jo = 1:Ng
        F(gi{jo},:) = squashme(F(gi{jo},:))+jo-1; % ff*(shmoo-min(shmoo(:)))/max(shmoo(:)) + jo; % also needs to span unit interval
    end

end
Ng = length(gi);
% set the colormap
if Ng <= 1,
    if strcmp(customcmap,'redyellowblue');
        customcmap = bengetcmap('redyellowblue',ng,0);
    else
        customcmap = gray(ng);
    end
    colormap(customcmap)
else
    cmap = colormap(bengetcmap('blues',ng,0,1));
    if Ng >= 2
        cmap = [cmap;bengetcmap('greens',ng,0,1)];
    end
    if Ng >= 3
        cmap = [cmap;bengetcmap('oranges',ng,0,1)];
    end
    if Ng >= 4
        cmap = [cmap;bengetcmap('purples',ng,0,1)];
    end
    if Ng >= 5
        cmap = [cmap;bengetcmap('reds',ng,0,1)];
    end
    if Ng >= 6
        cmap = [cmap;pink(ng)];
    end
    if Ng >= 7
        cmap = [cmap;gray(ng)];
    end
    if Ng >= 8
        cmap = [cmap;bengetcmap('yelloworangered',ng,0,1)];
    end
    if Ng >= 9
        cmap = [cmap;bengetcmap('purplebluegreen',ng,0,1)];
    end
    if Ng >= 10
        cmap = [cmap;bengetcmap('yellowgreenblue',ng,0,1)];
    end
    if Ng >= 11
        cmap = [cmap;bengetcmap('purpleblue',ng,0,1)];
    end
    if Ng >= 12
        cmap = [cmap;bengetcmap('purplered',ng,0,1)];
    end
    if Ng >= 13
        cmap = [cmap;bengetcmap('redpurple',ng,0,1)];
    end
    if Ng >= 14
        cmap = [cmap;bengetcmap('orangered',ng,0,1)];
    end
    if Ng >= 15
        cmap = [cmap;bengetcmap('yelloworangebrown',ng,0,1)];
    end
    if Ng >= 16
        cmap = [cmap;bengetcmap('greenblue',ng,0,1)];
    end
    if Ng >= 17
        cmap = [cmap;bengetcmap('bluepurple',ng,0,1)];
    end
    if Ng >= 18
        cmap = [cmap;bengetcmap('bluegreen',ng,0,1)];
    end
    if Ng >= 19
        disp('Unable to colour groups -- too many');
        cmap = bengetcmap('spectral',ng);
    end
    colormap(cmap)
end
% switch Ng
%     case {0,1}, colormap(bengetcmap('spectral',ng))
%     case 2, colormap([bengetcmap('redblue',ng);])
%     case 3, colormap([bengetcmap('redblue',ng);bengetcmap('purpleorange',ng);bengetcmap('browngreen',ng)])
%     case 4, colormap([bengetcmap('redblue',ng);bengetcmap('purpleorange',ng);bengetcmap('browngreen',ng)])
%     case 5, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng)])
%     case 6, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng)])
%     case 7, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng)])
%     case 8, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng);winter(ng)])
%     case 9, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng);winter(ng);copper(ng)])
%     otherwise, disp('too many groups to colour'), colormap(redbluecmap(ng));
% end
% switch Ng
%     case {0,1}, colormap(redbluecmap(ng))
%     case 2, colormap([gray(ng);pink(ng)])
%     case 3, colormap([gray(ng);pink(ng);spring(ng)])
%     case 4, colormap([gray(ng);pink(ng);spring(ng);summer(ng)])
%     case 5, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng)])
%     case 6, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng)])
%     case 7, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng)])
%     case 8, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng);winter(ng)])
%     case 9, colormap([gray(ng);pink(ng);spring(ng);summer(ng);bone(ng);hot(ng);cool(ng);winter(ng);copper(ng)])
%     otherwise, disp('too many groups to colour'), colormap(redbluecmap(ng));
% end


% Now reinvert F for operations, so that time series are again on the rows
if strcmp(metorts,'ops'); F = F'; end

% have to surround by zeros for an accurate pcolor:
pcolor([F zeros(size(F,1),1); zeros(1,size(F,2)+1)]);
shading flat

%%% Write the axes labels
%% fill axlabs with appropriately-formatted labels

% time series
if ~strcmp(labelm,'none')
    axlabs = BF_label_space(tsf,labelx);
    writeaxes(axlabs,'y')
end
    
% operations
if nm < 1000 % otherwise don't bother
    if ~strcmp(labelm,'none')
        axlabs = BF_label_space(mlab,labelx);
        writeaxes(axlabs,'x')
    end
end

%% Display the labels on the figure
set(gca,'FontSize',8)

end
