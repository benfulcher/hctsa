% TSQ_cluster
% 
% Reads in normalized data from HCTSA_N.mat, clusters the data matrix by
% reordering rows and columns, then saves the result as HCTSA_cl.mat
% 
% --INPUTS:
% ClusterMethRow: specifies the clustering method for rows/time series (default is 'linkage')
% ClusterParamsRow: specifies a cell of parameters specifying the clustering parameters,
%           including the distance metric, etc. (default is euclidean distances
%           and average linkage)
% ClusterMethCol: specifies the clustering method for columns/operations (default is the
%         same as rows)
% ClusterParamsCol: clustering settings for columns (default is correlation distances
%           and average linkage)
% SubSet: should be a cell with three components: the first should be either
%         'cl' or 'norm' for the row and column indicies provided in the 2nd and
%         3rd components of the cell. e.g., {'norm',[1,5,...],[]} will load in
%         HCTSA_N, and cluster the subset of rows [1,5,...] and all columns.
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

function TSQ_cluster(ClusterMethRow, ClusterParamsRow, ClusterMethCol, ClusterParamsCol, SubSet)

%% Check input arguments
if nargin < 1 || isempty(ClusterMethRow)
    ClusterMethRow = 'linkage';
end

% Clustering settings for rows
% (can specify the string 'none' for no clustering)
if nargin < 2
    ClusterParamsRow = {'euclidean','average',0,[],0}; % can change last 0 to 1 to save clustering info
end

% use same method as for rows
if nargin < 3 || isempty(ClusterMethCol)
    fprintf(1,'Clustering columns in the same way as for the rows: using %s\n',ClusterMethRow)
    ClusterMethCol = ClusterMethRow;
end

% clustering settings for columns
if nargin < 4 || isempty(ClusterParamsCol)
    ClusterParamsCol = {'corr','average',0,[],0}; % can change last 0 to 1 to save clustering info
%     ClusterParamsCol = ClusterParamsRow; % just specify the first one, and this will do the same
end

% Subsets -- only cluster a subset of the full data matrix
if nargin < 5
    fprintf(1,'No subsets, clustering the full object\n')
    SubSet = [];
end

if ~isempty(SubSet)
    if length(SubSet) ~= 3
        error('The subset should specify ''norm'' or ''cl'' and the subset.')
    elseif ~ismember(SubSet{1},{'norm','cl'})
        error('The first component of subset should be either ''norm'' or ''cl''.')
    end
end

%% Read in information from local files
if isempty(SubSet) || strcmp(SubSet{1},'norm')
    TheFile = 'HCTSA_N.mat';
else % subset of the clustered matrix
    TheFile = 'HCTSA_cl.mat';
end
wn = which(TheFile); % Check that HCTSA_N exists
if isempty(wn);
    error('%s not found.',TheFile);
end
fprintf(1,'Loading data from %s...',TheFile)
load('HCTSA_N.mat','TS_DataMat','TimeSeries','Operations','MasterOperations')
fprintf(1,' Done.\n')

% now all variables are by the 'cl' superscript names
% the data is in TS_DataMat

%% Implement subset behaviour
if ~isempty(SubSet)
    fprintf(1,'We are now implementing subset behaviour...\n')
    % subs is in the form {[rowrange],[columnrange]}; a cell of two vectors

    if ~isempty(SubSet{2}); % row subset
       r = SubSet{2};
       TS_DataMat = TS_DataMat(r,:);
       TimeSeries = TimeSeries(r);
       fprintf(1,'Reduced rows/timeseries from %u to %u according to the specified subset.\n',length(TimeSeries),length(r));
    end

    if ~isempty(SubSet{3}); % column subset
        r = SubSet{3};
        TS_DataMat = TS_DataMat(:,r);
        Operations = Operations(r);
        fprintf(1,'Reduced columns/operations from %u to %u according to the specified subset\n',length(Operations),length(r));
    end
end


%% Do the clustering
fprintf(1,'Clustering the full %u x %u data matrix.\n',length(TimeSeries),length(Operations))

% Cluster rows
if ~(ischar(ClusterMethCol) && ismember(ClusterMethRow,{'none','nothing'})) % can specify 'none' to do no clustering
    fprintf(1,'Clustering rows...\n'); tic
    [~, acgir] = TSQ_us_cluster(TS_DataMat,ClusterMethRow,ClusterParamsRow,'ts');
    fprintf(1,'Row clustering took %s.\n',BF_thetime(toc))
else
    acgir = {};
end

% Cluster columns
if ~(ischar(ClusterMethCol) && ismember(ClusterMethCol,{'none','nothing'})) && size(TS_DataMat,2)>1 % can specify 'none' to do no clustering
    fprintf(1,'Clustering columns...\n'); tic
    [~, acgic] = TSQ_us_cluster(TS_DataMat',ClusterMethCol,ClusterParamsCol,'mets');
    fprintf(1,'Column clustering took %s.\n',BF_thetime(toc))
else
    acgic = {};
end

%% Reorder output: TS_loc_cl is a reordering of F (TS_loc_N)
% get the permutation vectors ordr (reordering for rows) and ordc
% (reordering for columns)

if isempty(acgir)
    ordr = 1:size(TS_DataMat,1); % don't reorder at all
elseif iscell(acgir)
    ordr = vertcat(acgir{:});
else
    ordr = acgir;
end
if isempty(acgic);
    ordc = 1:size(TS_DataMat,2); % don't reorder at all
elseif iscell(acgic)
    ordc = vertcat(acgic{:});
else
    ordc = acgic;
end

% Reorder data matrix
TS_DataMat = TS_DataMat(ordr,ordc);

% Reorder time series metadata
TimeSeries = TimeSeries(ordr);

% Reorder operation metadata
Operations = Operations(ordc);

%% Save Output to file
% TS_loc_cl -- this is the clustered table with time series as rows and metrics as columns
fprintf(1,'Saving the clustered data as ''HCTSA_cl''...')
save('HCTSA_cl.mat','TS_DataMat','TimeSeries','Operations','MasterOperations')
fprintf(1,' Done.\n');

if length(ClusterParamsRow)>=5 && ClusterParamsRow{5}==1
    % reload the linkage saved from above and reorder it as per the new one
    disp('Reloading and reordering TS_guide_clinksr.mat!');
    load TS_guide_clinksr.mat R links
    R = squareform(R,'tomatrix');
    R = R(ordr,:); R = R(:,ordr);
    links12 = links(:,1:2);
    links12new = links12;
    for i = 1:length(R)
        links12new(links12==ordr(i)) = i;
    end
    links = [links12new,links(:,3)];
    R = squareform(R,'tovector');
    save TS_guide_clinksr.mat R links
end
if length(ClusterParamsCol)>=5 && ClusterParamsCol{5}==1
    % reload the linkage saved from above and reorder it as per the new one
    disp('Reloading and reordering TS_guide_clinksc.mat')
    load TS_guide_clinksc.mat R links
    R = squareform(R,'tomatrix');
    R = R(ordc,:); R = R(:,ordc);
    links12 = links(:,1:2);
    links12new = links12;
    for i = 1:length(R)
        links12new(links12==ordc(i)) =i;
    end
    links = [links12new,links(:,3)];
    R = squareform(R,'tovector');
    save TS_guide_clinksc.mat R links
end

% %% CRUDE: SAVE LINKAGE INFORMATION TO FILE
%     % In fact -- all we are saving are the reordered Rs calculated in the
%     % TSQ_us_cluster above. Do this using squareform(R), R = R(ord,:); R =
%     % R(:,ord)... But since we don't keep this in general, we can just do
%     % it this way -- just repeat the linkage with the reordered versions.
%     % May be a better idea to have a normalized and clustered version of
%     % the thing...
%     cparamsrr = cell(5,1); cparamscc = cell(5,1);
%     cparamsrr{5} = 1; % save to file
%     cparamscc{5} = 1; % save to file
%     for i = 1:4
%         if cparamsr>=i, cparamsrr{i} = cparamsr{i}; end
%         if cparamsc>=i, cparamscc{i} = cparamsc{i}; end
%     end
%     
%     % linkage and save to file
%     TSQ_us_cluster(TS_loc_cl,ClusterMethRow,cparamsrr,'ts');
%     TSQ_us_cluster(TS_loc_cl',cmethc,cparamscc,'mets');
%     
% %     % Do the linkage again on reordered versions: ... doesn't take too long, and helps with the indexing
% %     % Faster would be to just reorder based on ordc/ordr or something like that...
% %     % In fact, yes -- this is just the reordered Rs calculated in the
% %     % TSQ_us_cluster above. Do this using squareform(R), R = R(ord,:); R =
% %     % R(:,ord)... But since we don't keep this in general, we can just do
% %     % it this way.
% %     Rr = pdist(TS_loc_cl,dmth); linksr = linkage(Rr,lmth);
% %     Rc = pdist(TS_loc_cl',dmth); linksc = linkage(Rc,lmth);
%     
%     
%     
% %     % reorder
% %     [H,T,ord] = dendrogram(linksr,0); %,'colorthreshold','default'
% %     Rr = squareform(Rr); Rr = Rr(ord,:); Rr = Rr(:,ord);
% %     [H,T,ord] = dendrogram(linksc,0); %,'colorthreshold','default'
% %     Rc = squareform(Rc); Rc = Rc(ord,:); Rc = Rc(:,ord);
% %     Rr = squareform(Rr,'tovector'); % smaller for writing to file
% %     Rc = squareform(Rc,'tovector'); % smaller for writing to file
%     
%     % Save to file:
%     % TS_guide_clinks -- this contains the linkage information
% %     disp('Saving the linkage information as ''TS_guide_clinksr'' and ''TS_guide_clinksc''')
% %     save('TS_guide_clinksr.mat','Rr','linksr','-v7.3')
% %     save('TS_guide_clinksc.mat','Rc','linksc','-v7.3')
% end



%     function ord_cl = SUB_kmeanscl(F,k,distancemeasure,maxiterations)
%             a = kmeans; % specify the clustering model
%             a.k = k;
%             a.child = distance(distancemeasure); % set the distance measure
%             a.max_loops = maxiterations;
%             [r,a] = train(a,data(F)); % do the clustering
% 
%             % order first in terms of clusters
%             [rubbish ord_cl] = sort(r.X);
%             
%             % secondary: order by distance to cluster centre
%             cc = a.mu;
% %             ord_2 = zeros(size(F,1),1);
%             
%             for i=1:k
%                 ii = find(r.X==i);
%                 % in this subrange, order by distance to cluster centre
%                 dd = zeros(length(ii),1);
%                 for j=1:length(ii)
%                     dd(j) = calc(distance(distancemeasure),data(F(ii(j),:)),data(cc(i,:)));
%                     % evaluates the distance between each point and its
%                     % assigned cluster centre -- stores in vector dd
%                 end
%                 [rubbish ord_2] = sort(dd);
%                 ord_cl(ii) = ord_cl(ii(ord_2)); % reorder this segment by the secondary ordering
%             end
%     end


end