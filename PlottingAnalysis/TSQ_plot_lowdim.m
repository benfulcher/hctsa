% ------------------------------------------------------------------------------
% TSQ_plot_lowdim
% ------------------------------------------------------------------------------
%
% Calculates then plots a lower-dimensional feature-based representation of the
% data (e.g., using PCA).
% 
%----HISTORY:
% [Previously called TSQ_dimred]
% Ben Fulcher 31/3/2010 -- new classmeth option to specify classification
%                           method -- i.e., built in linear/quadratic; or
%                           svm-based method, etc.
% Ben Fulcher 18/4/2010 -- justus specifies whether to do PCA on the full matrix,
% 						   or just the groups specified in the given subset
% Ben Fulcher 13/7/2010 -- removed justus option!! Trying to clean up an
%                           inconsistency in labeling. Added TheData input
% Ben Fulcher 29/10/2010 -- added annotatep: can annotate time series
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

function TSQ_plot_lowdim(TheData,TsorOps,classmeth,showks,annotatep)

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 1 || isempty(TheData)
    TheData = 'norm';
    fprintf(1,'Getting data from HCTSA_N\n'); 
end

if nargin < 2 || isempty(TsorOps)
    fprintf(1,'Using for time series\n');
    TsorOps = 'ts';
end
if ~any(ismember(TsorOps,{'ops','ts'}));
    error('Specify either operations (''ops'') or time series (''ts'').');
end

% % Specify a keyword labeling of the data, kwgs:
% if nargin < 3
%     kwgs = {};
%     fprintf(1,'No assignment of the data? Ok up to you...\n');
% end
% 
% % Specify group indices, gi:
% if nargin < 4
%     gi = [];
%     fprintf(1,'Will obtain group indices from file\n');
% end

% if nargin < 5 || isempty(DimRedMethod)
%     DimRedMethod = 'pca';
%     fprintf(1,'Using PCA\n');
% end

if nargin < 3 || isempty(classmeth)
    classmeth = 'linclass';
    fprintf(1,'No discriminant\n'); 
end

if nargin < 4 || isempty(showks)
    showks = 1;
end

if nargin < 5 || isempty(annotatep)
    annotatep = struct('n',10);
end
if ~isstruct(annotatep)
    annotatep = struct('n',annotatep);
end

% ------------------------------------------------------------------------------
%% Load the data and group labeling from file
% ------------------------------------------------------------------------------
if strcmp(TheData,'cl') || strcmp(TheData,'norm') 
    % Retrive data from local files
    switch TheData
    case 'cl'
        TheDataFile = 'HCTSA_cl.mat';
    case 'norm'
        TheDataFile = 'HCTSA_N.mat';
    end
    fprintf(1,'Loading data and grouping information from %s...',TheDataFile);
    load(TheDataFile,'TS_DataMat');
    if strcmp(TsorOps,'ts')
        load(TheDataFile,'Operations')
        DimensionLabels = {Operations.Name}; clear Operations % We just need their names
        if isstruct(annotatep) || length(annotatep) > 1 || annotatep > 0
            load(TheDataFile,'TimeSeries')
            DataLabels = {TimeSeries.FileName};
            data_ids = [TimeSeries.ID];
            TimeSeriesData = {TimeSeries.Data};
            if isfield(TimeSeries,'Group')
                load(TheDataFile,'GroupNames')
                DataGroups = [TimeSeries.Group];
            else
                fprintf(1,'\n');
                error('No groups assigned -- Use TSQ_LabelGroups.')
            end
            clear('TimeSeries'); % we no longer need you
        end
    else
        load(TheDataFile,'TimeSeries')
        DimensionLabels = {TimeSeries.FileName}; clear TimeSeries
    end
    fprintf(1,' Loaded.\n');
else
    % The user provided data yourself
    TS_DataMat = TheData;
    if isfield(annotatep,'olab')
        DimensionLabels = annotatep.olab;
    else
        DimensionLabels = {};
    end
end

if strcmp(TsorOps,'ops')
    % Take the transpose of the input data matrix for operations
    TS_DataMat = TS_DataMat';
end


GroupIndices = BF_ToGroup(DataGroups);
% if isempty(gi)
%     gi = SUB_autolabelQ(kwgs,TsorOps,TheData);
% end
% CheckEmpty = cellfun(@isempty,gi);
% if any(CheckEmpty)
%     error('No keywords found for: %s . Exiting.',kwgs{find(CheckEmpty,1)})
% end
% if (size(kwgs,2) == 2) && (size(kwgs,1) > 1); % specified subsets of each keyword
%    kwgs = kwgs(:,1);
% end
NumGroups = length(GroupIndices); % Number of groups

% ------------------------------------------------------------------------------
%% Do the dimensionality reduction
% ------------------------------------------------------------------------------

% Matlab's build-in PCA
% Sort it so that when choose different set of keywords the output is consistent
% There's a strange thing in princomp that can give different scores when the
% input rows are in a different order. The geometry is the same, but they're reflected
% relative to different orderings
fprintf(1,'Calculating principal components of the %u x %u data matrix...', ...
                    size(TS_DataMat,1),size(TS_DataMat,2));
[pc,score,latent] = princomp(TS_DataMat);
fprintf(1,' Done.\n');
PercVar = round(latent/sum(latent)*1000)/10; % Percentage of variance explained (1 d.p.)

% Work out the main contributions to each principle component
featlabel = cell(2,2); % Feature label :: proportion of total [columns are PCs)
tolabel = cell(2,1);
ngcontr = min(length(DimensionLabels),2);  % the 2 greatest contributions to the first principle component
for i = 1:2
    [s1, ix1] = sort(abs(pc(:,i)),'descend');
    featlabel{i,1} = DimensionLabels(ix1(1:ngcontr));
    featlabel{i,2} = round(abs(s1(1:ngcontr)/sum(abs(s1)))*100)/100;
    for j = 1:ngcontr
        tolabel{i} = sprintf('%s%s (%f), ',tolabel{i},featlabel{i,1}{j},featlabel{i,2}(j));;
    end
    tolabel{i} = tolabel{i}(1:end-2);
end

% ------------------------------------------------------------------------------
% Plot this two-dimensional representation of the data
% ------------------------------------------------------------------------------

% Move all 2d plotting to TSQ_plot_2d
% extras.TS_DataMat = score;

NameString = 'PC';
for i = 1:2
    DataInfo.labels{i} = sprintf('%s%u (%f%%) : %s',NameString,i,PercVar(i),tolabel{i});
end

% function TSQ_plot_2d(Features,DataInfo,TrainTest,annotatep,keepksdensities,lossmeth,extras)

DataInfo.GroupNames = GroupNames;
DataInfo.GroupIndices = GroupIndices;
DataInfo.DataLabels = DataLabels;
DataInfo.TimeSeriesData = TimeSeriesData;

TSQ_plot_2d(score(:,1:2),DataInfo,{},annotatep,showks,classmeth);

% if ischar(TheData)
% else
%     % Hopefully tsf is specified
%     extras.tsf = annotatep.tsf;
%     TSQ_plot_2d(score(:,1:2),DataInfo,{},annotatep,showks,classmeth,extras);
% end


% % set up sc => first 2 principal components grouped
% gi_r = gi;
% sc = cell(NumGroups,1);
% for i = 1:NumGroups
%     sc{i} = score(gi_r{i},1:2); % include first two principle components
% end

% % define colours for plotting
% if NumGroups>8
%     c = bengetcmap('set3',NumGroups,1);
% else
%     c = bengetcmap('set1',NumGroups,1);
% end
% 
% figure('color','w')
% 
% %% Plot distributions
% 
% if showks
%     % top one
%     subplot(4,4,1:3); hold on; box('on')
%     for i=1:NumGroups
%         plot_ks(sc{i}(:,1),c{i},0)
%     end
%     set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[])
% 
%     % right side one
%     subplot(4,4,[8 12 16]); hold on; box('on')
%     for i=1:NumGroups
%         plot_ks(sc{i}(:,2),c{i},1)
%     end
%     set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[])
% end
% 
% 
% %% Plot 1,2 Principal components
% if showks
%     subplot(4,4,[5:7,9:11,13:15]);
% end    
% box('on'); hold on;
% switch classmeth
%     case {'linear','quadratic','svm','none'}
%         for i = 1:NumGroups
%             % points
%             plot(sc{i}(:,1),sc{i}(:,2),'.','MarkerSize',8,'color',c{i})
%         end
%         for i = 1:NumGroups
%             % cluster centres
%             cc = mean(sc{i}); %median(sc{i});
%             plot(cc(1),cc(2),'o','color',c{i},'MarkerFaceColor',c{i}+(1-c{i})*0.5,...
%                             'MarkerSize',10,'LineWidth',2);
%         end
% %         if shadenice
% %             % partition the space into regions and colour based on this
% %             nbins = 25; % number of partitions along each axis
% %             c = cell(4,1); % for four groups
% %             c{1} = bengetcmap('blues',5,1);
% %             c{2} = bengetcmap('greens',5,1);
% %             c{3} = bengetcmap('oranges',5,1);
% %             c{4} = bengetcmap('purples',5,1);
% %             
% %             % make 2d histogram for all groups involved
% %             hist2s = zeros(NumGroups,nbins,nbins);
% %             for i = 1:NumGroups
% %                 v1 = sc{i}(:,1);
% %                 v2 = sc{i}(:,2);
% %                 edgesi = givemeedges('range',v1,nbins);
% %                 [ni bini] = histc(v1, edgesi);
% %                 edgesj = givemeedges('range',v2,nbins);
% %                 [nj binj] = histc(v2, edgesj);
% % 
% %                 % CREATE JOINT HISTOGRAM
% %                 % we have the edges in each dimension: edgesi, and edgesj
% %                 hist2s{i} = zeros(nbins);
% %                 for i2 = 1:nbins
% %                     for j2 = 1:nbins
% %                         hist2s(i,i2,j2) = sum(bini==i2 & binj==j2);
% %                     end
% %                 end
% %             end
% %             
% %             % Now find winner in each box
% %             winners = zeros(nbins);
% %             for i = 1:nbins
% %                 for j = 1:nbins
% %                     winners(i,j) = find(max(hist2s(:,i,j)));
% %                 end
% %             end
% %             
% %             % partition the space
% %             totalhits = squeeze(sum(hist2s,1));
% %             totalhits = totalhits(:);
% %             thresh = quantile(totalhits,linspace(0,1,5));
% %             
% %             % now plot
% %             for i = 1:nbins
% %                 for j = 1:nbins
% %                     
% %                 end
% %             end
% %         end
%             
%     case 'gmm'
% %         Nclasses = 2;
%         Nclasses = NumGroups;
%         X = vertcat(sc{:}); % instead of score
%         i1 = 1:length(sc{1}); i2 = length(sc{1})+1:size(X,1);
%         options = statset('Display','final');
%         gm = gmdistribution.fit(X,Nclasses,'Replicates',5,'Options',options);
% %         scatter(sc{1}(:,1),sc{1}(:,2),12,P(gi_r{1},1),'+')
% %         plot(gm.mu(1,1),gm.mu(1,2),'ok'); % the mean of this gaussian mixture
% %         hold on
% %         scatter(sc{2}(:,1),sc{2}(:,2),12,P(gi_r{2},1),'o')
% %         plot(gm.mu(2,1),gm.mu(2,2),'or'); % the mean of this gaussian mixture
% %         hold off
% 
%         hold on
%         
% %         % first the data:
%         P = posterior(gm,X);
% %         ai1 = find(P(:,1)<=0.5); ai2 = find(P(:,1)>0.5); % autoassignments based on posterior
%         idx = cluster(gm,X); % get cluster assignments
%         ai1 = find(idx==1); ai2 = find(idx==2);
%         
% %         % Option 1: plot automatic compared to labelled assignments:
% %         % plot correctly assigned
% %         C1 = intersect(i1,ai1); C2 = intersect(i2,ai2);
% %         I1 = intersect(i1,ai2); I2 = intersect(i2,ai1);
% %         if length(I1)>length(C1),
% %             cc1 = C1; C1 = I1; I1 = cc1;
% %             cc2 = C2; C2 = I2; I2 = cc2;
% %             c1 = 'b'; c2 = 'r';
% %         else
% %             c1 = 'r'; c2 = 'b';
% %         end
% %         scatter(X(C1,1),X(C1,2),9,'r','+') 
% %         scatter(X(C2,1),X(C2,2),9,'b','x')
% %         
% %         % plot incorrectly assigned
% %         scatter(X(I1,1),X(I1,2),20,'r','o')
% %         scatter(X(I2,1),X(I2,2),20,'b','o')
% 
%         % Option 2: just plot automatic assignments by posterior
% %         scatter(X(ai1,1),X(ai1,2),10,P(ai1,1),'+')
% %         scatter(X(ai2,1),X(ai2,2),10,P(ai2,1),'x')
% %         c1 = [1,0,1]; c2 = [0,1,1]; % use these colours
%         
%         % Option 3: plot automatic assignments in red/blue
%         scatter(X(ai1,1),X(ai1,2),10,'r','+')
%         scatter(X(ai2,1),X(ai2,2),10,'b','x')
%         c1 = 'r'; c2 = 'b';
%         
% %         keyboard
% %         scatter(X(i1,1),X(i1,2),30,'r','.') 
% %         scatter(X(i2,1),X(i2,2),30,'b','.')
% %         
% %         % now the cluster centres (means of gaussian mixture components)
% 
% 
%         plot(gm.mu(1,1),gm.mu(1,2),'o','color','k','MarkerFaceColor',c1,'MarkerSize',12,'LineWidth',2);
%         plot(gm.mu(2,1),gm.mu(2,2),'o','color','k','MarkerFaceColor',c2,'MarkerSize',12,'LineWidth',2);
%         
%         % now the contours
%         xr = get(gca,'xlim');yr = get(gca,'ylim');
% %         xr = [-1,1];yr=[-1,1];
%         ezcontour(@(x,y)pdf(gm,[x y]),xr,yr,400);
%         colormap('cool')
%         
% end
% 
% %% -- plot a classify boundary
% if NumGroups == 2
%     switch classmeth
%         case {'linear','quadratic'}
%             modeorder = classmeth; % or quadratic -- or add in an svm classifier?
% 
%             xlim=get(gca,'XLim'); ylim=get(gca,'YLim');
%             group = givemegroup(gi_r);
% %             group=zeros(size(score,1),1);
% %             for i=1:NumGroups, group(gi_r{i})=i; end
%             [X,Y] = meshgrid(linspace(xlim(1),xlim(2),200),linspace(ylim(1),ylim(2),200));
%             X = X(:); Y = Y(:);
%             [C,err,P,logp,coeff] = classify([X Y],[score(:,1) score(:,2)],group, modeorder);
% 
%             hold on;
%             %     gscatter(X,Y,C,'rb','o',1,'off');
%             K = coeff(1,2).const; L = coeff(1,2).linear;
%             if strcmp(modeorder,'linear')
%                 Q = zeros(2,2);
%             else
%                 Q = coeff(1,2).quadratic;
%             end
%             f = sprintf('0 = %g+%g*x+%g*y+%g*x^2+%g*x.*y+%g*y.^2',K,L,Q(1,1),Q(1,2)+Q(2,1),Q(2,2));
%             h2 = ezplot(f,[xlim(1) xlim(2) ylim(1) ylim(2)]);
%             set(h2,'LineStyle','--','color','k','LineWidth',3)
%         case 'svm'
%             disp('SVM BOUNDARY?! HELP!')
%             keyboard
%             d = makeitdata(score(:,1:2),gi_r);
%             [r,a] = train(svm,d);
%             plot(a{1})
%     end
% end
% 
% 
% %% Labelling
% switch DimRedMethod
%     case 'pca'
%         xlabel(['PC1 (' num2str(PercVar(1)) '%) : ' tolabel{1}],'interpreter','none')
%         ylabel(['PC2 (' num2str(PercVar(2)) '%) : ' tolabel{2}],'interpreter','none')
%         title('');
%     case 'nmf'
%         xlabel(['NMF1: ' tolabel{1}],'interpreter','none')
%         ylabel(['NMF2: ' tolabel{2}],'interpreter','none')
%         title(['Mean reconstruction error: ' num2str(meanreconstructionerror)]);
% end
% 
% % if showks
% %     ylabel('');
% % end
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
% 
% 
% % Set up legend
% legs = cell(NumGroups,1);
% for i=1:NumGroups, legs{i}=[kwgs{i} ' (' num2str(length(gi{i})) ')']; end
% legend(legs,'interpreter','none');
% 
% % quote loss in title
% if NumGroups>1
%     try
%         loss = TSQ_givemeloss(kwgs,gi,{'knn',{1}},{'cv',{5,1}},'class_loss',F);
%         title(['knn(1) loss = ' num2str(round(loss*100)) ' %   with ' num2str(size(F,2)) ' features'])
%     catch
%         disp('No Spider')
%     end
% end
% 
% %% Which original axes (features) align most closely to the new principal
% %% components?
% % keyboard
% % biplot(pc(:,1:2),'Scores',score(:,1:2),'VarLabels',...
% % 		olab')
% 
% %% ANNOTATE
% if isstruct(annotatep)
%     if isfield(annotatep,'maxL')
%         maxL = annotatep.maxL;
%     else
%         maxL = 300; % length of annotated time series segments
%     end
%     if isfield(annotatep,'uinput')
%         uinput = annotatep.uinput;
%     else
%         uinput = 0; % user input points rather than randomly chosen
%     end
%     if isfield(annotatep,'fdim')
%         fdim = annotatep.fdim;
%     else
%         fdim = []; % set default below
%     end
%     annotatep = annotatep.n;
% else % specify parameters as vector
%     if length(annotatep)>=2
%         maxL = annotatep(2);
%     else
%         maxL = 300; % length of annotated time series segments
%     end
% 
%     if length(annotatep)>=3
%         uinput = annotatep(3);
%     else
%         uinput = 0; % user input points rather than randomly chosen
%     end
%     annotatep = annotatep(1);
%     fdim = [];
% end
% 
% 
% if annotaten>0 % define xy: coordinates of points in the space
%     xy = cell(NumGroups,1);
%     for i = 1:NumGroups
%         xy{i} = [sc{i}(:,1),sc{i}(:,2)];
%     end
% end
% 
% pwidth = diff(get(gca,'xlim')); % plot width
% pheight = diff(get(gca,'ylim')); % plot height
% if isempty(fdim), fdim = [0.23,0.07]; end % width, height
% lwidth = 0.7; % line width for time series
% textann = 'tsid'; % no annotations (also 'filename','tsid')
% plotcircle = 1; % magenta circle around annotated points
% 
% alreadypicked = zeros(annotaten,2); % record those already picked
% if ~uinput % points to annotate are randomly picked
%     if annotaten == length(tsf) % annotate all
%         disp('annotate all')
%        for j = 1:annotaten
%            thegroup = find(cellfun(@(x)ismember(j,x),gi));
%            alreadypicked(j,1) = thegroup;
%            alreadypicked(j,2) = find(gi{thegroup}==j);
%        end
%     else
%     %     thegroups = randi(NumGroups,[annotaten,1]); % random groups
%         alreadypicked(:,1) = round(linspace(1,NumGroups,annotaten));
%         randperms = cellfun(@(x)randperm(length(x)),gi,'UniformOutput',0);
%         counters = ones(NumGroups,1);
%         for j=1:annotaten
%             alreadypicked(j,2) = randperms{alreadypicked(j,1)}(counters(alreadypicked(j,1))); % random element of the group
%             counters(alreadypicked(j,1)) = counters(alreadypicked(j,1))+1;
%         end
%     end
% end
% for j = 1:annotaten
% % %     user input:
%     if uinput
%         point = ginput(1);
%         iplot = benginpclosest(xy,point); % find closest actual point to input point
%         thegroup = iplot(1); % want this group
%         itsme = iplot(2); % and this index
%         alreadypicked(j,:) = [thegroup,itsme];
%     else
%         thegroup = alreadypicked(j,1);
%         itsme = alreadypicked(j,2);
%     end
%     
% %     keyboard
%     if j>1 && any(sum(abs(alreadypicked(1:j-1,:)-repmat(alreadypicked(j,:),j-1,1)),2)==0) % same already been picked
%         continue; % don't plot this again
%     end
%     
%     plotpoint = xy{thegroup}(itsme,:);
%     fn = tsf{gi{thegroup}(itsme)}; % filename of timeseries to plot
%     ts = dlmread(fn);
%     ts = ts(1:min(maxL,end));
%     if plotcircle
%         plot(plotpoint(1),plotpoint(2),'om'); % plot magenta circle around target point
%     end
%     if strcmp(textann,'filename')
%         % annotate text with filename:
%         text(plotpoint(1),plotpoint(2)-0.01*pheight,fn,'interpreter','none','FontSize',8);
%     elseif strcmp(textann,'tsid');
%         % annotate text with ts_id:
%         text(plotpoint(1),plotpoint(2)-0.01*pheight,num2str(ts_ids_keep(gi{thegroup}(itsme))),'interpreter','none','FontSize',8);
%     end
%     plot(plotpoint(1)-fdim(1)*pwidth/2+linspace(0,fdim(1)*pwidth,length(ts)),...
%             plotpoint(2)+fdim(2)*pheight*(ts-min(ts))/(max(ts)-min(ts)),...
%             '-','color',c{thegroup},'LineWidth',lwidth);
% end
% 
% 
% % alreadypicked = zeros(annotaten,2); % record those already picked
% % for j = 1:annotaten
% % % %     user input:
% % %     point = ginput(1);
% % %     iplot = benginpclosest(xy,point); % find closest actual point to input point
% % %     thegroup = iplot(1); % want this group
% % %     itsme = iplot(2); % and this index
% %     % or random point:
% %     thegroup = randi(NumGroups); % random group
% %     itsme = randi(length(gi{thegroup})); % random element of the group
% %     
% %     alreadypicked(j,:) = [thegroup,itsme];
% % %     keyboard
% %     if j>1 && any(sum(alreadypicked(1:j-1,:)-repmat(alreadypicked(j,:),j-1,1),2)==0) % same already been picked
% %         continue; % don't plot this again
% %     end
% %     
% %     plotpoint = xy{thegroup}(itsme,:);
% %     fn = tsf{gi{thegroup}(itsme)}; % filename of timeseries to plot
% %     ts = dlmread(fn);
% %     ts = ts(1:min(maxL,length(ts)));
% %     pwidth = diff(get(gca,'xlim')); % plot width
% %     pheight = diff(get(gca,'ylim')); % plot height
% % %     plot(plotpoint(1),plotpoint(2),'om'); % plot magenta circle around
% % %     the target point
% %     % annotate text with filename:
% % %     text(plotpoint(1),plotpoint(2)-0.01*pheight,fn,'interpreter','none','FontSize',8);
% %     % annotate text with ts_id:
% % %     text(plotpoint(1),plotpoint(2)-0.01*pheight,num2str(ts_ids_keep(gi{thegroup}(itsme))),'interpreter','none','FontSize',8);
% %     plot(plotpoint(1)+linspace(0,0.20*pwidth,length(ts)),...
% %         plotpoint(2)+0.05*pheight*(ts-min(ts))/(max(ts)-min(ts)),...
% %         '-','color',c{thegroup},'LineWidth',1);
% % end
% 
% 
% function plot_ks(v,c,swap)
%     % vector v is the vector of a given group
%     % c is the color
%     [f x]=ksdensity(v,linspace(min(v),max(v),1000),'function','pdf');
%     r=zeros(length(v),1);
%     for m=1:length(v); r(m)=find(x>=v(m),1,'first'); end
%     r=sort(r,'ascend'); r=r(round(linspace(1,length(r),50))); % crop to just 50 points across the range
%     if swap
%         plot(f,x,'color',c);
%         plot(f(r),x(r),'.','color',c)
%     else
%         plot(x,f,'color',c);
%         plot(x(r),f(r),'.','color',c)
%     end
% end

end