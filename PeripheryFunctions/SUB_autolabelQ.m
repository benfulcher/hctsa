function gi = SUB_autolabelQ(kwgs,metorts,norcl,kworlab,subset)
% INPUTS:
% (*) The keyword groups, a cell of strings: kwgs
% (*) Whether keywords are for metrics or time series: metorts='mets' or 'ts'
% (*) Whether to retrive from normal or clustered store: norcl='orig', 'norm', or 'cl'
% (*) [opt] kworlab: can specify 'lab' for the actual specific names

% OUTPUTS:
% (*) gi: the indicies corresponding to each keyword in kwgs

% Ben Fulcher August 2009
%       New options: can specify kwgs with numbers: e.g.,
%       kwgs={'space',100;'',200;'medical',0;...} [0 means all]
%       [blank '' label means anything: selected at random from all time series]
% ------------------------------------------------------------------------
% Ben Fulcher 24/3/2010 -- changed to SUB_autolabelQ from SUB_autolabel2 to
%                           use the SQL server to retrieve.
% Easier to get the keywords matching those specified from the SQL Server
% and then looking for overlaps with those present
%
% Ben Fulcher 30/9/2010 -- added subset option to specify a custom subset
% to label (in that index system)

%% (0) Check inputs
if nargin < 1 || isempty(kwgs), disp('you must specify labels'); return; end
if nargin < 2 || isempty(metorts), disp('Assuming time series'); metorts = 'ts'; end
if nargin < 3 || isempty(norcl), disp('Assuming you want to retrieve from clustered index system'); norcl = 'cl'; end
if nargin < 4 || isempty(kworlab), kworlab = 'kw'; end % look up keywords by default
if nargin < 5, subset = []; end % don't subset by default

if ~ismember(metorts,{'mets','ts'})
    disp('SUB_autolabel: Specify ''mets'' or ''ts''. Exiting.'), return
end
if ~isstruct(norcl) && ~ismember(norcl,{'orig','norm','cl'})
    disp('SUB_autolabel: Specify ''orig'', ''norm'', or ''cl''. Exiting.'), return
end
if ~ismember(kworlab,{'kw','lab'})
    disp('SUB_autolabel: Specify ''kw'' or ''lab''. Exiting.'), return
end
if ischar(kwgs); kwgs={kwgs,0}; end


%% (1) The keywords are loaded from guides and stored in kws
if isstruct(norcl)
    kws = norcl.kws;
    idsO = norcl.idsO;
else
    if strcmp(metorts,'ts'),
        switch norcl
            case 'orig'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides.mat tskw ts_ids_keep
                    kws = tskw; idsO = ts_ids_keep;
                else
                    load TS_loc_guides.mat tsf ts_ids_keep
                    kws = tsf; idsO = ts_ids_keep;
                end
            case 'norm'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides_N.mat tskwn ts_ids_keepn
                    kws = tskwn; idsO = ts_ids_keepn;
                else
                    load TS_loc_guides_N.mat tsfn ts_ids_keepn
                    kws = tsfn; idsO = ts_ids_keepn;
                end
            case 'cl'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides_cl.mat tskwcl ts_ids_keepcl
                    kws = tskwcl; idsO = ts_ids_keepcl;
                else
                    load TS_loc_guides_cl.mat tsfcl ts_ids_keepcl
                    kws = tsfcl; idsO = ts_ids_keepcl;
                end
        end
    else % operations
        switch norcl
            case 'orig'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides.mat mkw m_ids_keep
                    kws = mkw; idsO = m_ids_keep;
                else
                    load TS_loc_guides.mat mlab m_ids_keep
                    kws = mlab; idsO = m_ids_keep;
                end
            case 'norm'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides_N.mat mkwn m_ids_keepn
                    kws = mkwn; idsO = m_ids_keepn;
                else
                    load TS_loc_guides_N.mat mlabn m_ids_keepn
                    kws = mlabn; idsO = m_ids_keepn;
                end
            case 'cl'
                if strcmp(kworlab,'kw')
                    load TS_loc_guides_cl.mat mkwcl m_ids_keepcl
                    kws = mkwcl; idsO = m_ids_keepcl;
                else
                    load TS_loc_guides_cl.mat mlabcl m_ids_keepcl
                    kws = mlabcl; idsO = m_ids_keepcl;
                end
        end
    end
end

if ~isempty(subset)
    kws = kws(subset);
    idsO = idsO(subset);
end


%% Option (1)
% % Get all from database and return intersection with ids
% 
% % Let's get the matching indicies from the database
% % Open database -> dbc
% dbc = SQL_opendatabase(dbname);
% 
% % Get matching ids from the database for each keyword
% % Use TSQ_getids
% Ng = size(kwgs,1);
% gi = cell(Ng,1); % stores the indicies (in the clustering index system) in a vector
% 
% % This is one option: get all matching IDS, and find intersection. Another
% % option is to look through the obtained keywords from the local guide -- 
% % i.e., do it all in MATLAB. This was originally done, but now the keywords
% % are comma-delimited rather than cells of cells.
% for i = 1:Ng
%     allmatchingids = TSQ_getids(metorts,[],kwgs(i,:),{},[],[],dbname);
%     gi{i} = arrayfun(@(x)find(idsO==x),intersect(allmatchingids, idsO));
%     if isempty(gi{i})
%         disp(['Nothing found for ' kwgs{i,1} ' :-(']);
%     end
% end
% 
% % Done -- close database
% mysql_dbexecute(dbc, 'closeall');


%% Option 2
% Search in SQL query for intersection

%% Option 3
% Convert to cell and to it all within MATLAB
if ~all(cellfun(@ischar,kwgs(:))) % have specified numbers of each
    kwnums = horzcat(kwgs{:,2}); % just the number of each part
    kwgs = kwgs(:,1)'; % just the keyword parts, a cell of strings
else
    kwnums=zeros(length(kwgs),1); % include all of each keyword
end
Ng = length(kwgs); % number of groups
% ndata = length(kws); % number of time series/operations
kws = SUB_cell2cellcell(kws);

if strcmp(kworlab,'kw') % look for groups of keywords
%     wally = zeros(ndata,1);
    % ; corresponding to the keywords in kwgs
%     disp('Automatic labeling from the guide initiated.')
    for jo=1:Ng
        if ~isempty(kwgs{jo}) % collect time series with this keyword
%             for i=1:ndata, wally(i)=any(ismember(kws{i},kwgs{jo})); end
            wally = cellfun(@(x)any(ismember(kwgs{jo},x)),kws);
            gi{jo} = find(wally);
            if isempty(gi{jo})
                disp(['Problem with group ' kwgs{jo} ': no matches found.']);
%                 return
            end
            if kwnums(jo)~=0 && kwnums(jo)<length(gi{jo}) % take a random subset
                rperm = randperm(length(gi{jo}));
                gi{jo} = gi{jo}(rperm(1:kwnums(jo)));
            end
        else % take a certain number of random time series
             % integer: retrieve this many: in randomorder
            rperm = randperm(length(kws));
            kwgs{jo} = 'OTHERS';
            gi{jo} = [];
            notkwgs = kwgs(setxor(1:Ng,jo));
            for i=1:length(kws)
                if all(~ismember(notkwgs,kws{rperm(i)}))
                    gi{jo} = [gi{jo}; rperm(i)];
                    if length(gi{jo})==kwnums(jo)
                        break
                    end
                end
            end
        end
    end
    disp('Autolabelling complete. Now that wasn''t so bad, was it?!');
    

else % find the indicies of specific labels (no ability for numbers of each;
      % -- there should be only one
    gi=zeros(Ng,1);
    for jo=1:Ng
        gi(jo) = strmatch(kwgs{jo},kws,'exact');
        if isempty(gi(jo))
            disp(['Problem with ' kwgs{jo} ': no matches found.']);
%             return
        end
    end
end

end

