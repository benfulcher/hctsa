function TSQ_prepared(ts_ids_keep, m_ids_keep, getwhat, dbname, brawninputs, doinone)

% This function prepares files for subsequent calculation/visualization methods.
% It takes as input a set of constraints on the time series and metrics to include
% and output the relevant subsection of the time series-metric table (taken from 
% STORE_DAT) as TS_loc and also a matlab file containints the time series labels, lengths,
% keywords, and preprocessings (this is TS_guide_ts.mat)
% and the same for metric function names, labels, post-processings, and master/pointer labels
% which is stored as TS_guide_met.mat
% It also initializes the STORE_data matrix with missing values along rows and columns corresponding
% to new time series/metrics being added for the first time to the store
% 
% HISTORY:
% 29/10/2009: added changed functionality; only saves STORE backup if index system has changed.
% 12/11/2009: added logging functionality; logs to the local directory a summary of what's been done.
% 9/1/2010: added trimming functionality; trim to some number based on more than just random (and do this within mySQL syntax)
% 10/1/2010: add masterpull option; by default pulls all other metrics which masters call and adds to 
% 				the operation list (for calculation purposes this is sensible)
% 11/1/2010: delegates the subsetting role to TSQ_getids -- should call this either in argument (e.g., 
% 				TSQ_prepared(TSQ_getids('ts',...),TSQ_getids('mets',...))) or call them first and call these vectors
% 				to TSQ_prepared
% 12/1/2010: tried to make so only retrieve bits required for calculation if forcalc flag is enabled. If forcalc = [],
% 				 will retrieve everything. If zero will retrieve all empty entries. Otherwise will retrieve up to that
% 				 number of entries. (at this stage, QualityCode=NULL; eventually could add fatal error ones)
% 12/5/2010: added doinone -- database will either make many connections (suitable for small no. of rows), or do it in one
% 			 big connection (strain on java heap space for large
% 			 retrievals)

%% INPUTS:
% 1) lenr: contrain included time series by length [minimum_length maximum_length] (1x2 vector)
% 2) kyesc: constrain included time series by keyword -- nx2 cell
% 			(i) keyword
% 			(ii) how many of that keyword (0==all available)
% 3) kno: keywords NOT to include (cell of strings)
% 4) myesc: metrics to include (nx2 cell)
% 			(i) keyword
% 			(ii) how many of that keyword (0==all available)
% 5) mno: metrics not to include (cell of strings)
% 7) masterpull [opt]: default 1 -- stops
% doinone -- can specify to 1 to retrieve all on a single database connection

%% OUTPUTS (to file):
% 1) TS_loc.mat, which contains the matrix TS_loc: the portion of the storage file that
% 				   matches the supplied constaints.
% 2) TS_loc_guides.mat, which contains information about the time series in TS_loc:
%				 () tsf: filenames of the files containing the time series data (cell of strings)
%				 () tsl: the length of each time series (vector)
%				 () tskw: the keywords of each time series (cell of string cells)
%				 () tsprep: the type of detrending required (vector of integer labels) [[NOW REDUNDANT]]
%				 () tsmap: the mapping of the time series in the portion TS_loc to STORE_data (integer vector)
%				 () nts: the number of time series (scalar integer)
% 				 () mf: the matlab code to call for each metric (cell of strings)
% 				 			(or the relevant part of a master structure output if a pointer)
% 				 () mlab: labels of the metrics; different to their calling function (cell of strings)
% 				 () mpostp: post-processings specific to each metric (vector of integer labels)
% 				 () mkw: keywords for each metric (cell of string cells)
% 				 () mtyp: the type of metric: i.e., single (S), or pointer (P)
% 				 () mmap: mapping from the metrics in TS_loc to their position in STORE_data (integer vector)
% 				 () nm: the number of metrics (scalar integer)
% 				 () mMf: the matlab code to call (cell of strings)
% 				 () mMl: master metric labels, for pointers to point to (cell of strings)
% 				 () nmM: the number of master metrics (scalar integer)


%%% FOREPLAY
%% Check inputs -- set defaults
if nargin < 2;
	error('You must provide at least two inputs!');
end
if nargin < 3
	getwhat = 'all'; % retrieve full sets of things, not just the empty entries in the database
end
if nargin < 4
	dbname = []; % Use default database
end
if nargin < 5 || isempty(brawninputs)
	brawninputs = [1, 1]; % log and parallelize -- only relevant if getwhat isn't empty
end
if nargin < 6 || isempty(doinone)
	doinone = 0; % make seperate connections so as not to overwhelm java heap space
end

%% METHOD 1: Get entries of results table for local MATLAB matrix (1 row/column at a time)
% we have ts_ids_keep and m_ids_keep
% To put in a matrix with rows (time series) and columns (metrics)
% Could do one big query and then reform to a matrix, but I'll do it row-by-row
% In fact this is faster for some reason than doing a big query (method 2)

% make sure ts_ids_keep and m_ids_keep are column vectors
if size(ts_ids_keep,2) > size(ts_ids_keep,1), ts_ids_keep = ts_ids_keep'; end
if size(m_ids_keep,2) > size(m_ids_keep,1), m_ids_keep = m_ids_keep'; end
ts_ids_keep = sort(ts_ids_keep,'ascend'); % sorted by ascending id
m_ids_keep = sort(m_ids_keep,'ascend'); % sorted by ascending id
ts_ids_keep_string = bencat(ts_ids_keep,',');
m_ids_keep_string = bencat(m_ids_keep,',');
nts = length(ts_ids_keep);
nm = length(m_ids_keep);

if nts == 0 || nm == 0
	fprintf(1,'Oops. It seems there''s nothing to do! Either no time series or no operations\n');
	return
end

% Intialize matrices
% Start as Infs to distinguish unwritten entries after database pull
TS_loc = ones(nts,nm)*Inf; % outputs
TS_loc_ct = ones(nts,nm)*Inf; % calculation times
TS_loc_q = ones(nts,nm)*Inf; % output quality label

% Open database connection
[dbc,dbname] = SQL_opendatabase(dbname);

% Tell me about it
fprintf(1,'We have %u time series and %u operations to retrieve from %s\n',nts,nm,dbname);
fprintf(1,'Filling and saving to local matricies TS_loc, TS_loc_ct, TS_loc_q from Results table in %s\n',dbname);

bundlesize = 5; % retrieve information about this many time series per database query

switch getwhat
case 'null'
    fprintf(1,'Retrieving NULL elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
case 'all' 
    fprintf(1,'Retrieving all elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
case 'error'
    fprintf(1,'Retrieving error elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
end

bundles = (1:bundlesize:length(ts_ids_keep));
nits = length(bundles); % number of iterations of the loop
times = zeros(nits,1); % record time for each iteration
% First select the data
for i = 1:nits
    ittic = tic; % start timing the iteration

    ii = (bundles(i):1:min(bundles(i)+bundlesize-1,length(ts_ids_keep))); % indicies for this iteration
    ts_ids_now = ts_ids_keep(ii); % a range of ts_ids to retrieve in this iteration
    basestring = sprintf(['SELECT ts_id, m_id, Output, CalculationTime, QualityCode FROM Results WHERE ' ...
                    	'ts_id IN (%s) AND m_id IN (%s)'],bencat(ts_ids_now,','),m_ids_keep_string);
    switch getwhat
    case 'null'
    	SelectString = [basestring, ' AND QualityCode IS NULL'];
    case 'error'
    	SelectString = [basestring, ' AND QualityCode = 1'];
    case 'all'
        SelectString = basestring;
    end
	[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString); % retrieve the bundlesize from the database
    
    % Check results look ok:
    if ~isempty(emsg)
        fprintf(1,'Error retrieving outputs from %s???\n',dbname);
        fprintf(1,'%s\n',emsg)
        keyboard
    end
    if size(qrc) == 0
        fprintf(1,'No data to retrieve for ts_id = (%s)\n',bencat(ts_ids_now,','));
    end
    
	% Convert empty entries to NaNs
	qrc(cellfun(@isempty,qrc)) = {NaN};
    
    % Fill rows of local matrix
    if strcmp(getwhat,'all') % easy in this case
        TS_loc(ii,:) = reshape(vertcat(qrc{:,3}),length(ii),nm);
        TS_loc_ct(ii,:) = reshape(vertcat(qrc{:,4}),length(ii),nm);
        TS_loc_q(ii,:) = reshape(vertcat(qrc{:,5}),length(ii),nm);
    else
        % we need to match retrieved indicies with local indicies
        ix = arrayfun(@(x)find(ts_ids_now==x,1),vertcat(qrc{:,1}));
        iy = arrayfun(@(x)find(m_ids_keep==x,1),vertcat(qrc{:,2}));
        for k = 1:length(ix) % fill it one entry at a time
            TS_loc(ix(k),iy(k)) = qrc{k,3};
            TS_loc_ct(ix(k),iy(k)) = qrc{k,4};
            TS_loc_q(ix(k),iy(k)) = qrc{k,4};
        end
    end
    
    % Note time taken for this iteration
	times(i) = toc(ittic);
	if mod(i,floor(nits/10)) == 0 % tell the user 10 times
		fprintf(1,'Approximately %s remaining...\n',benrighttime(mean(times(1:i))*(nts-i)));
	end
end
fprintf(1,'Local files filled from %s in %u iterations. Took %s altogether.\n',dbname,nits,benrighttime(sum(times)));
	
if ismember(getwhat,{'null','error'})    
    % We only want to keep rows and columns with (e.g., NaNs) in them...
    switch getwhat
    case 'null'
        keepme = isnan(TS_loc); % NULLs in database
    	fprintf(1,'Filtering so that data only contains rows/columns containing at least one NULL entry\n');
    case 'error'
        keepme = (TS_loc_q == 1); % error codes in database
    	fprintf(1,'Filtering so that data only contains rows/columns containing at least one error entry\n');
    end
    
	% Time Series
    keepi = (sum(keepme,2) > 0); % there is at least one entry to calculate in this row
    if sum(keepi)==0
    	fprintf(1,'After filtering, there are no time series remaining! Exiting...\n'); return
	elseif sum(keepi) < nts
		fprintf(1,'Cutting down from %u to %u time series\n',nts,sum(keepi));
		ts_ids_keep = ts_ids_keep(keepi); nts = length(ts_ids_keep);
		ts_ids_keep_string = bencat(ts_ids_keep,',');
		TS_loc = TS_loc(keepi,:); TS_loc_ct = TS_loc_ct(keepi,:); TS_loc_q = TS_loc_q(keepi,:); % update local data stores
	end
	
	% Operations
    keepi = (sum(keepme,1) > 0); % there is at least one entry to calculate in this column
	if sum(keepi) == 0
    	fprintf(1,'After filtering, there are no operations remaining! Exiting...\n'); return
    elseif sum(keepi) < nm
		fprintf(1,'Cutting down from %u to %u operations\n',nm,sum(keepi));
		m_ids_keep = m_ids_keep(keepi); nm = length(m_ids_keep);
		m_ids_keep_string = bencat(m_ids_keep,',');
		TS_loc = TS_loc(:,keepi); TS_loc_ct = TS_loc_ct(:,keepi); TS_loc_q = TS_loc_q(:,keepi);
	end    
end

%% Fill out guides
% Retrieve Time Series Metadata
SelectString = sprintf('SELECT FileName, Keywords, Length FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_keep_string);
[tsinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
tsf = tsinfo(:,1); % timeseries filenames
tskw = tsinfo(:,2); % timeseries keywords
tsl = vertcat(tsinfo{:,3}); % timeseries lengths

% Retrieve Operations Metadata
SelectString = sprintf('SELECT Code, OpName, Keywords FROM Operations WHERE m_id IN (%s)',m_ids_keep_string);
[minfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
mcode = minfo(:,1); % operation filenames
mlab = minfo(:,2); % operation labels
mkw = minfo(:,3); % operation keywords
% mpoint = vertcat(minfo{:,4}); % pointer? -- [redundant if you have mlink (retrieved later)]

% Now get master info
% (i) Which masters are implicated?
SelectString = ['SELECT mop_id, MasterLabel, MasterCode FROM MasterOperations WHERE mop_id IN ' ...
				'(SELECT DISTINCT mop_id FROM MasterPointerRelate WHERE m_id IN (' m_ids_keep_string '))'];
[masterinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
if ~isempty(emsg)
    fprintf(1,'Error retrieving Master information...\n'); keyboard
else
    if ~isempty(masterinfo) % there are masters in out midst
        Mmid = vertcat(masterinfo{:,1});
        Mmlab = masterinfo(:,2);
        Mmcode = masterinfo(:,3);
    else
        % No master functions implicated in this subset
        Mmid = [];
        Mmlab = {};
        Mmcode = {};
    end
end

% (ii) Get master link information
SelectString = ['SELECT mop_id FROM (SELECT m_id FROM Operations WHERE m_id IN (' m_ids_keep_string ')) AS T1 ' ...
					'LEFT JOIN MasterPointerRelate ON T1.m_id = MasterPointerRelate.m_id'];
[masterlink,~,~,emsg] = mysql_dbquery(dbc,SelectString);
empties = find(cellfun(@isempty,masterlink));
if ~isempty(empties), masterlink(empties) = {0}; end
mlink = vertcat(masterlink{:});

% (iii) renumber mlink to refer to elements of Mm rather than the indicies in Master index system
nM = length(Mmid); % number of master metrics
rch = cell(nM,1);
for i = 1:nM
    rch{i} = (mlink==Mmid(i));
end
for i = 1:nM
    mlink(rch{i}) = i; % rename with index rather than the mop_id
end

% Close database connection
SQL_closedatabase(dbc)

fprintf(1,'Saving local versions of the data...:\n');
save('TS_loc.mat','TS_loc','-v7.3')
fprintf(1,'TS_loc (data)')

save('TS_loc_ct.mat','TS_loc_ct','-v7.3')
fprintf(1,', TS_loc_ct (calculation times)')
	
save('TS_loc_q.mat','TS_loc_q','-v7.3')
fprintf(1,', TS_loc_q (quality codes)')

save('TS_loc_guides.mat','m_ids_keep','ts_ids_keep','tsf','tskw','tsl','mcode','mlab','mkw','mlink','Mmid','Mmlab','Mmcode','-v7.3');
fprintf(1,', TS_loc_guides (information guides).\n')

% Display how many entries need to be calculated
tocalculate = sum(isnan(TS_loc_q(:)) | TS_loc_q(:)==1);
fprintf(1,'There are %u entries (= %5.2f%%) to calculate in TS_loc (%ux%u)\n',tocalculate,tocalculate/nm/nts*100,nts,nm);

% reply = input(['Shall I move on to TSQ_brawn now to calculate and then write back? [yes, ''n'' to stop]'],'s');
% if ~isempty(forcalc)
%     fprintf(1,'Automatically moving on to TSQ_brawn now to calculate the missing entries and then write them back to the database\n');
%     TSQ_brawn(brawninputs(1),brawninputs(2),dbname);
% end

end