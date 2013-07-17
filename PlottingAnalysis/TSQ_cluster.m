function TSQ_cluster(cmethr, cparamsr, cmethc, cparamsc, subs) % ,actuallycluster
%%% TS_cluster %%%
% Takes in T_loc_N (the post-processed data), clusters the object, saving the result as
% TS_loc_cl (the data) and TS_guide_cl_ts (time series labels) and TS_guide_cl_met (metric labels)
% actuallycluster is an optional boolean; if 0, will skip the clustering steps;
% and just do the filtering of bad columns/rows bit.
% optional input howstrict, if 1, will remove any columns with NaNs

% (*)subs(*): should be a cell with three components: the first should be either
% 'cl' or 'norm' for the row and column indicies provided in the 2nd and
% 3rd components of the cell. e.g., {'norm',[1,5,...],[]} will load in
% TS_loc_N, and cluster the subset of rows [1,5,...] and all columns.
% 
% Ben Fulcher 9/12/2009 -- Updated to use mySQL storage system
% Ben Fulcher 1/4/2010 -- Rehauled to take advantage of spider clustering
%                           routines
% Ben Fulcher 7/6/2010 -- added seperate clustering options for rows and
%                           columns
% Ben Fulcher 15/6/2010 -- added savelinkage option -- if doing linkage
%                           clustering, whether to save the linkage
%                           information to file. A bit hacky to append this
%                           like this. Removed on 17/6/2010.

%% Check input arguments

if nargin < 1 || isempty(cmethr)
    cmethr = 'linkage';
%     cparams = {'euclidean','single',0}; % {dmth,lmth,showdend}
%   these defaults are set within the method itself
end

% clustering settings for rows
% (can specify the string 'none' for no clustering)
if nargin < 2
    cparamsr = {'euclidean','average',0,[],0}; % can change last 0 to 1 to save clustering info
end

% use same method as for rows
if nargin < 3 || isempty(cmethc)
    fprintf(1,'Clustering in the same way as for the columns: using %s\n',cmethr)
    cmethc = cmethr;
    cparamsc = cparamsr; % just specify the first one, and this will do the same
end

% clustering settings for columns
if nargin < 4 || isempty(cparamsc)
    cparamsc = {'corr','average',0,[],0}; % can change last 0 to 1 to save clustering info
%     cparamsc = cparamsr; % just specify the first one, and this will do the same
end

% subsets -- only cluster a subset of the full data matrix
if nargin < 5
    fprintf(1,'No subsets, clustering the full object\n')
    subs = [];
end

if ~isempty(subs)
    if length(subs) ~= 3
        error('The subset should specify ''norm'' or ''cl'' and the subset.')
    elseif ~ismember(subs{1},{'norm','cl'})
        error('The first component of subset should be either ''norm'' or ''cl''.')
    end
end

%% Read in information from local files
fprintf(1,'Loading in local files...')

if isempty(subs) || strcmp(subs{1},'norm')
    fprintf(1,' TS_loc_N, ...')
    wn = which('TS_loc_N.mat'); % check it exists
    if isempty(wn); fprintf(1,'\n'), error('TS_loc_N not found -- please run TSQ_normalize'); end
    load TS_loc_N.mat TS_loc_N % this is the normalized local data file -- we cluster on the normalized values
    F = TS_loc_N;
    clear TS_loc_N
    
    load TS_loc_guides_N.mat m_ids_keepn ts_ids_keepn tsfn tskwn tsln mcoden mlabn mkwn mlinkn Mmid Mmlab Mmcode
    % change names to clustered versions with superscript 'cl'
    % m_ids_keep, ts_ids_keep, tsf, tskw, tsl, mcode, mlab, mkw, mlink, Mmid, Mmlab, Mmcode

    % () time series ()
    ts_ids_keepcl = ts_ids_keepn; tsfcl = tsfn; tskwcl = tskwn; tslcl = tsln;

    % () metrics ()
	m_ids_keepcl = m_ids_keepn; mcodecl = mcoden; mlabcl = mlabn; mkwcl = mkwn; mlinkcl = mlinkn; % mpointcl = mpointn;
else % subset of the clustered matrix
    load TS_loc_cl.mat TS_loc_cl
    load TS_loc_guides_cl.mat % this is the relevant information about the time series and operations
    F = TS_loc_cl;
    clear TS_loc_cl
end
fprintf(1,' loaded.\n')

% now all variables are by the 'cl' superscript names
% the data is in F

%% Implement subset behaviour
if ~isempty(subs)
%     if strcmp(subs{1},'norm') % SUBSET OF TS_loc_N
        disp('We are now implementing subset behaviour')
        % subs is in the form {[rowrange],[columnrange]}; a cell of two vectors
        
        if ~isempty(subs{2}); % row subset
           r = subs{2};
           disp(['Subsetting rows/timeseries: from ' num2str(length(tsfcl)) ' to ' num2str(length(r))]);
           F = F(r,:);
           ts_ids_keepcl = ts_ids_keepcl(r); tsfcl=tsfcl(r); tskwcl=tskwcl(r); tslcl=tslcl(r);
        end
        
        if ~isempty(subs{3}); % column subset
            r = subs{3};
            disp(['Subsetting columns/metrics: from ' num2str(length(mkwcl)) ' to ' num2str(length(r))]);
            F = F(:,r);
            m_ids_keepcl = m_ids_keepcl(r); mcodecl = mcodecl(r); mlabcl=mlabcl(r); mkwcl=mkwcl(r); mlinkcl = mlinkcl(r); % mpointcl = mpointcl(r); 
        end
else
    disp('Preparing to cluster the full (TS_loc_N) matrix')
end


%% Do the clustering
	
% Cluster rows
if ~(ischar(cmethc) && ismember(cmethr,{'none','nothing'})) % can specify 'none' to do no clustering
    disp('clustering rows...'); tic
    [~, acgir] = TSQ_us_cluster(F,cmethr,cparamsr,'ts');
    disp(['row clustering took ' BF_thetime(toc)])
else
    acgir = {};
end

% Cluster columns
if ~(ischar(cmethc) && ismember(cmethc,{'none','nothing'})) && size(F,2)>1 % can specify 'none' to do no clustering
    disp('clustering columns...'); tic
    [~, acgic] = TSQ_us_cluster(F',cmethc,cparamsc,'mets');
    disp(['column clustering took ' BF_thetime(toc)])
else
    acgic = {};
end


%% Reorder output: TS_loc_cl is a reordering of F (TS_loc_N)
% get the permutation vectors ordr (reordering for rows) and ordc
% (reordering for columns)

if isempty(acgir)
    ordr = 1:size(F,1); % don't reorder at all
elseif iscell(acgir)
    ordr = vertcat(acgir{:});
else
    ordr = acgir;
end
if isempty(acgic);
    ordc = 1:size(F,2); % don't reorder at all
elseif iscell(acgic)
    ordc = vertcat(acgic{:});
else
    ordc = acgic;
end

% reorder data matrix
TS_loc_cl = F(ordr,ordc);

% Reorder row guides
ts_ids_keepcl = ts_ids_keepcl(ordr); tsfcl=tsfcl(ordr); tskwcl=tskwcl(ordr); tslcl=tslcl(ordr);

% Reorder column guides
m_ids_keepcl = m_ids_keepcl(ordc); mcodecl = mcodecl(ordc); mlabcl=mlabcl(ordc); mkwcl=mkwcl(ordc);
mlinkcl = mlinkcl(ordc); % mpointcl = mpointcl(ordc); 



%% Save Output to file
% TS_loc_cl -- this is the clustered table with time series as rows and metrics as columns
disp('Saving the clustered data as ''TS_loc_cl''')
save('TS_loc_cl.mat','TS_loc_cl')

% TS_guide_cl_ts, TS_guide_cl_met -- this contains all the clustered time series 
% and metric information, respectively.
% Note that the clustering is only done on 'good' metrics and so nmcl<=nm
disp('Saving guides: ''TS_loc_guides_cl.mat''...');
save('TS_loc_guides_cl.mat','m_ids_keepcl','ts_ids_keepcl','tsfcl','tskwcl','tslcl','mcodecl','mlabcl','mkwcl','mlinkcl','Mmid','Mmlab','Mmcode','-v7.3');


if length(cparamsr)>=5 && cparamsr{5}==1
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
if length(cparamsc)>=5 && cparamsc{5}==1
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
%     TSQ_us_cluster(TS_loc_cl,cmethr,cparamsrr,'ts');
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