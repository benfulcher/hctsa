function TSQ_prepared(ts_ids_keep, m_ids_keep, tolog, forcalc, dbname, brawninputs, doinone)

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
% 				 number of entries. (at this stage, Quality=NULL; eventually could add fatal error ones)
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
% 6) tolog [opt]: whether or not to log to a local file a summary of the input/output. (default = 0 -- don't log)
% 7) masterpull [opt]: default 1 -- stops
% doinone -- can specify to 1 to retrieve all on a single database connection

%% OUTPUTS (to file):
% 1) TS_loc.mat, which contains the matrix TS_loc: the portion of the storage file that
% 				   matches the supplied constaints.
% 2) TS_guide_ts.mat, which contains information about the time series in TS_loc:
%				 (i) tsf: filenames of the files containing the time series data (cell of strings)
%				 (ii) tsl: the length of each time series (vector)
%				 (iii) tskw: the keywords of each time series (cell of string cells)
%				 (iv) tsprep: the type of detrending required (vector of integer labels) [[NOW REDUNDANT]]
%				 (v) tsmap: the mapping of the time series in the portion TS_loc to STORE_data (integer vector)
%				 (vi) nts: the number of time series (scalar integer)
% 3) TS_guide_met.mat, which contains information about the metrics in TS_loc:
% 				 (i) mf: the matlab code to call for each metric (cell of strings)
% 				 			(or the relevant part of a master structure output if a pointer)
% 				 (ii) mlab: labels of the metrics; different to their calling function (cell of strings)
% 				 (iii) mpostp: post-processings specific to each metric (vector of integer labels)
% 				 (iv) mkw: keywords for each metric (cell of string cells)
% 				 (v) mtyp: the type of metric: i.e., single (S), or pointer (P)
% 				 (vi) mmap: mapping from the metrics in TS_loc to their position in STORE_data (integer vector)
% 				 (vii) nm: the number of metrics (scalar integer)
% 				 (viii) mMf: the matlab code to call (cell of strings)
% 				 (ix) mMl: master metric labels, for pointers to point to (cell of strings)
% 				 (x) nmM: the number of master metrics (scalar integer)


%%% FOREPLAY
%% Check inputs -- set defaults
if nargin<2;
	disp('You must provide at least 2 inputs! And no, I''m not asking nicely!!');
	return
end
if nargin<3 || isempty(tolog)
	disp('Not logging TSQ_prepared');
	tolog = 0;  % no log created
end
if nargin<4
	forcalc = []; % retrieve full sets of things, not just the empty entries in the database
end
if nargin<5
	dbname = []; % Use default
end
if nargin<6 || isempty(brawninputs)
	brawninputs = [1 1]; % log and parallelize -- only relevant if forcalc isn't empty
end
if nargin<7 || isempty(doinone)
	doinone = 0; % make seperate connections so as not to overwhelm java heap space
end

% if nargin<7 || isempty(masterpull)
% 	masterpull = 1;
% 	disp(['Pulling in pointers by default']);
% end

% if isempty(lenr)
% 	lenr=[200,30000];
% 	disp(['Setting default length constraints: ' num2str(lenr(1)) '--' num2str(lenr(2))])
% end

% chind=0; % changed index system boolean.
%          % If stays 0 throughout, no need to resave or save backup
%          % If 1, have added a time series or metric and will need to
%          % resave the storage files


% dbc = SQL_opendatabase(dbname); % Open SQL connection

% length constraint in words
% ss_tsl = ['TimeSeries.Length BETWEEN ' num2str(lenr(1)) ' AND ' num2str(lenr(2))];

%% DISPLAY WHAT'S BEEN CHOSEN
% selectstring = ['SELECT FileName FROM TimeSeries WHERE ts_id IN (' ...
% 					bencat(ts_ids_keep) ')'];
% [ts_filenames,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% ts_filenames
% 
% 
% selectstring = ['SELECT OpName FROM Operations WHERE m_id IN (' ...
% 					bencat(m_ids_keep) ')'];
% [m_names,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% m_names
% 
% keyboard

%% METHOD 1: Get entries of results table for local MATLAB matrix (1 row/column at a time)
% we have ts_ids_keep and m_ids_keep
% To put in a matrix with rows (time series) and columns (metrics)
% Could do one big query and then reform to a matrix, but I'll do it row-by-row
% In fact this is faster for some reason than doing a big query (method 2)

ts_ids_keep = sort(ts_ids_keep,'ascend'); % sorted by ascending id
m_ids_keep = sort(m_ids_keep,'ascend'); % sorted by ascending id
ts_ids_keep_string = bencat(ts_ids_keep,',');
m_ids_keep_string = bencat(m_ids_keep,',');
nts = length(ts_ids_keep);
nm = length(m_ids_keep);

% Tell me about it
disp(['We have ' num2str(nts) ' time series and ' num2str(nm) ' operations to retrieve']);

if nts==0 || nm==0
	disp('Nothing to do! Enough of this.');
	return
end

disp('Filling local matricies TS_loc, TS_loc_ct, TS_loc_q from Results table in database');

TS_loc = zeros(nts,nm); % outputs
TS_loc_ct = zeros(nts,nm); % calculation times
TS_loc_q = zeros(nts,nm); % output quality label

if ~isempty(forcalc)
	disp('Retrieving elements to calculate from the database. Please be patient... I''m timing.');
	tic
	% Just retrieve bits that need to be calculated
	% May be faster to break it up and do each ts_id seperately... (this way the amount of data in each transfer is limited)
	if forcalc==0
		selectstring = ['SELECT ts_id, m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ')' ...
							' AND m_id IN (' m_ids_keep_string ' ) AND Quality IS NULL'];
	else
		selectstring = ['SELECT ts_id, m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ')' ...
							' AND m_id IN (' m_ids_keep_string ' ) AND Quality IS NULL LIMIT ' num2str(forcalc)];
	end
	dbc = SQL_opendatabase(dbname,0); % Open SQL connection silently
	[qrc,~,~,emsg] = mysql_dbquery(dbc,selectstring);
	SQL_closedatabase(dbc); % close connections
	if ~isempty(emsg)
		disp('Error retrieving outputs from database -- Exiting'); keyboard
	else
		ngot = size(qrc,1);
		if ngot==0
			disp('NOTHING TO RETRIEVE?!!');
			SQL_closedatabase(dbc)
			return
		else
			disp(['Retrieved ' num2str(ngot) ' elements!!! Took ' benrighttime(toc)]);
		end
	end
	
	% Turn empty entries into NaNs.
	qrc(cellfun(@isempty,qrc)) = {NaN};
	
	% segment outputs; give them sensible names
	Gts_ids = vertcat(qrc{:,1});
	Gm_ids = vertcat(qrc{:,2});
	Goutput = horzcat(qrc{:,3});
	Gcalctime = horzcat(qrc{:,4});
	Gquality = horzcat(qrc{:,5});
	
	% augment Gid_pairs into indicies
	Gindex_ts = Gts_ids;
	for i=1:nts
		Gindex_ts(Gts_ids == ts_ids_keep(i)) = i; % switch from the ts_id to the index
	end
	Gindex_m = Gm_ids;
	for i=1:nm
		Gindex_m(Gm_ids == m_ids_keep(i)) = i; % switch from the m_id to the index
	end
	% Gindex_pairs = [Gindex_ts,Gindex_m]
	
	% We have the index pairs, we now must just write these into the matricies
	% First, initiate matricies with Infs (means bogus -- unattained)
	
	disp('Filling indicies...');
	tic
	TS_loc(:)=Inf; TS_loc_ct(:)=Inf; TS_loc_q = Inf;
	for i=1:ngot
		TS_loc(Gindex_ts(i),Gindex_m(i)) = Goutput(i);
		TS_loc_ct(Gindex_ts(i),Gindex_m(i)) = Gcalctime(i);
		TS_loc_q(Gindex_ts(i),Gindex_m(i)) = Gquality(i);
	end
	disp(['Local stores filled. Took ' benrighttime(toc)]);
	
	
	% ()()() Filter down TS_loc_ and ts_ids_keep, m_ids_keep ()()()
	disp('Filtering unused time series and operations');
	tic
	% Time Series
	Utsids = unique(Gindex_ts); % rows with things to calculate in them (will be rows of Infs in TS_loc_)
	nUtsids = length(Utsids);
	if nUtsids < nts
		disp(['Cutting down from ' num2str(nts) ' to ' num2str(nUtsids) ' time series with calculation potential']);
		ts_ids_keep = ts_ids_keep(Utsids);
		TS_loc = TS_loc(Utsids,:);
		TS_loc_ct = TS_loc_ct(Utsids,:);
		TS_loc_q = TS_loc_q(Utsids,:);
	end
	
	% Operations
	Umids = unique(Gindex_m); % rows with things to calculate in them (will be rows of Infs in TS_loc_)
	nUmids = length(Umids);
	if nUmids < nm
		disp(['Cutting down from ' num2str(nm) ' to ' num2str(nUmids) ' operations with calculation potential']);
		m_ids_keep = m_ids_keep(Umids);
		TS_loc = TS_loc(:,Umids);
		TS_loc_ct = TS_loc_ct(:,Umids);
		TS_loc_q = TS_loc_q(:,Umids);
	end
	
	if nUtsids < nts || nUmids < nm % there has been cutting done. Update things.
		% ts_ids_keep = sort(ts_ids_keep,'ascend'); % sorted by ascending id
		% m_ids_keep = sort(m_ids_keep,'ascend'); % sorted by ascending id
		ts_ids_keep_string = bencat(ts_ids_keep,',');
		m_ids_keep_string = bencat(m_ids_keep,',');
		nts = length(ts_ids_keep);
		nm = length(m_ids_keep);
	end
	
	disp(['Filtering complete: took ' benrighttime(toc)]);
	
	% keyboard
	
else % retrieve everything in the range given (not just empty bits)
	% choose whether to do by rows or columns
	% dorows = 1;
	% if dorows
	
	
	% Fill a row at a time
	% Opening and closing individual connections helps to avoid a massive java heap space error... Hopefully...
	times = zeros(nts,1);
	if doinone
		dbc = SQL_opendatabase(dbname);
	end
	for i=1:nts
		tic
		selectstring = ['SELECT m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i)) ...
							' AND m_id IN (' m_ids_keep_string ' )'];
		if ~doinone, dbc = SQL_opendatabase(dbname,0); end % Open SQL connection silently
		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
		if ~doinone, SQL_closedatabase(dbc); end % close connection
		if ~isempty(emsg)
			disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
		end
	
		% Convert empty entries to NaNs
		qrc(cellfun(@isempty,qrc)) = {NaN};
		% empties = find(cellfun(@isempty,qrc));
		% if ~isempty(empties), qrc(empties) = {NaN}; end
	
		% give sensible names
		Gm_id = horzcat(qrc{:,1});
		Goutput = horzcat(qrc{:,2});
		Gcalctime = horzcat(qrc{:,3});
		Gquality = horzcat(qrc{:,4});
	
		% make sure m_ids are in ascending order to match
		% (should be, if table is ordered properly...)
		[S_m_ids ix] = sort(Gm_id);
	
		% Check matches
		if ~all(S_m_ids' - m_ids_keep == 0)
			disp(['Problem with ' num2str(ts_ids_keep(i))]);
			keyboard
		end
	
		Goutput = Goutput(ix);
		Gcalctime = Gcalctime(ix);
		Gquality = Gquality(ix);

		TS_loc(i,:) = Goutput;
		TS_loc_ct(i,:) = Gcalctime;
		TS_loc_q(i,:) = Gquality;

		times(i) = toc;
		if mod(i,floor(nts/10))==0
			disp(['Approximately ' benrighttime(mean(times(1:i))*(nts-i)) ' remaining!']);
		end
	end
	if doinone
		SQL_closedatabase(dbc);
	end

	disp(['Took ' benrighttime(sum(times)) ' altogether']);

	% else % fill a column at a time
	% 	times = zeros(nm,1);
	% 	for i=1:nm
	% 		tic
	% 		selectstring = ['SELECT Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ' ) ' ...
	% 						'AND m_id = ' num2str(m_ids_keep(i))];
	% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
	% 		if ~isempty(emsg)
	% 			disp(['Error retrieving outputs from STORE at ' num2str(m_ids_keep(i)) ' -- Exiting']); keyboard
	% 		end
	% 	
	% 		% Convert empty entries to NaNs
	% 		qrc(cellfun(@isempty,qrc)) = {NaN};
	% 		% empties = find(cellfun(@isempty,qrc));
	% 		% if ~isempty(empties)
	% 		% 	qrc(empties) = {NaN};
	% 		% end
	% 
	% 		TS_loc(:,i) = vertcat(qrc{:,1});
	% 		TS_loc_ct(:,i) = vertcat(qrc{:,2});
	% 		TS_loc_q(:,i) = vertcat(qrc{:,3});
	% 
	% 		times(i) = toc;
	% 
	% 		if mod(i,floor(nm/10))==0
	% 			disp(['Approximately ' benrighttime(mean(times(1:i))*(nm-i)) ' remaining!']);
	% 		end
	% 	end
end


%% METHOD 1B: GET NON-NULL ENTRIES A ROW/COLUMN AT A TIME...
% % We have ts_ids_keep and m_ids_keep to put in a matrix with rows (time series) and columns (metrics)
% % Could do one big query and then reform to a matrix, but I'll do it row-by-row
% % In fact this is faster for some reason than doing a big query (method 2)
% disp('Filling local matricies TS_loc, TS_loc_ct, TS_loc_q from Results table in database');
% ts_ids_keep = sort(ts_ids_keep,'ascend'); % make sure in ascending order
% m_ids_keep = sort(m_ids_keep,'ascend'); % make sure in ascending order
% ts_ids_keep_string = bencat(ts_ids_keep,',');
% m_ids_keep_string = bencat(m_ids_keep,',');
% 
% TS_loc = zeros(nts,nm); % outputs
% TS_loc_ct = zeros(nts,nm); % calculation times
% TS_loc_q = zeros(nts,nm); % output quality label
% 
% % choose whether to do by rows or columns
% dorows = 1;
% 
% if dorows
% 	times = zeros(nts,1);
% 	for i=1:nts
% 		tic
% 		selectstring = ['SELECT m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i)) ...
% 							' AND m_id IN (' m_ids_keep_string ' ) AND Output IS NULL'];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% 		if ~isempty(emsg)
% 			disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
% 		end
% 		
% 		if ~isempty(qrc)
% 			% Convert empty entries to NaNs
% 			empties = find(cellfun(@isempty,qrc));
% 			if ~isempty(empties), qrc(empties) = {NaN}; end
% 		
% 			% give sensible names
% 			Gm_id = horzcat(qrc{:,1});
% 			Goutput = horzcat(qrc{:,2});
% 			Gcalctime = horzcat(qrc{:,3});
% 			Gquality = horzcat(qrc{:,4});
% 		
% 			% make sure m_ids are in ascending order to match
% 			% (should be, if table is ordered properly...)
% 			[S_m_ids ix] = sort(Gm_id);
% 			Goutput = Goutput(ix);
% 			Gcalctime = Gcalctime(ix);
% 			Gquality = Gquality(ix);
% 		
% 			% Put into required spots; since only retrieved a subset
% 			mapper = zeros(length(S_m_ids),1); % maps from sorted m_ids from query (S_m_ids) to entries in local store
% 			for j=1:length(S_m_ids)
% 				mapper(j) = find(m_ids_keep == S_m_ids(j),1,'first');
% 			end
% 		
% 			% fill entries :: put nonzero bits in
% 			% (***) Using this method, the quality will be artificially zero for all non-null outputs
% 			% (***) Using this method, the calculation time will be artificially zero for all non-null outputs
% 			% (***) Using this method, fatal errors stored with output 0 will not be retrieved and subsequently recalculated
% 			% (***) Using this method, real outputs, properly-calculated output will be artificially zero.
% 			TS_loc(i,mapper) = Goutput;
% 			TS_loc_ct(i,mapper) = Gcalctime;
% 			TS_loc_q(i,mapper) = Gquality;
% 		% else
% 		% 	disp(['No entries need retrieving for ' num2str(ts_ids_keep(i))]);
% 		end
% 
% 		times(i) = toc;
% 	
% 		if mod(i,floor(nts/10))==0
% 			disp(['Approximately ' benrighttime(mean(times(1:i))*(nts-i)) ' remaining!']);
% 		end
% 	end
% 
% 	disp(['Took ' benrighttime(sum(times)) ' altogether']);
% 	
% else % fill a column at a time
% 	times = zeros(nm,1);
% 	for i=1:nm
% 		tic
% 		selectstring = ['SELECT Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ' ) ' ...
% 						'AND m_id = ' num2str(m_ids_keep(i))];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% 		if ~isempty(emsg)
% 			disp(['Error retrieving outputs from STORE at ' num2str(m_ids_keep(i)) ' -- Exiting']); keyboard
% 		end
% 		
% 		% Convert empty entries to NaNs
% 		empties = find(cellfun(@isempty,qrc));
% 		if ~isempty(empties)
% 			qrc(empties) = {NaN};
% 		end
% 	
% 		TS_loc(:,i) = vertcat(qrc{:,1});
% 		TS_loc_ct(:,i) = vertcat(qrc{:,2});
% 		TS_loc_q(:,i) = vertcat(qrc{:,3});
% 
% 		times(i) = toc;
% 	
% 		if mod(i,floor(nm/10))==0
% 			disp(['Approximately ' benrighttime(mean(times(1:i))*(nm-i)) ' remaining!']);
% 		end
% 	end
% 	
% 	
% end



%% METHOD 2: Get entries of results table for local MATLAB matrix
% % we have ts_ids_keep and m_ids_keep
% % To put in a matrix with rows (time series) and columns (metrics)
% % Could do one big query and then reform to a matrix, but I'll do it row-by-row
% disp('Filling local matricies TS_loc, TS_loc_ct, TS_loc_q');
% nts = length(ts_ids_keep); % number of time series
% disp([num2str(nts) ' timeseries']);
% nm = length(m_ids_keep); % number of metrics
% disp([num2str(nm) ' operations']);
% ts_ids_keep_string = bencat(ts_ids_keep,',');
% m_ids_keep_string = bencat(m_ids_keep,',');
% 
% TS_loc = zeros(nts,nm); % outputs
% TS_loc_ct = zeros(nts,nm); % calculation times
% TS_loc_q = zeros(nts,nm); % output quality label
% 
% % Retieve all data in one go from database
% disp(['Retrieving data from database...']); tic
% selectstring = ['SELECT ts_id, m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ') ' ...
% 					'AND m_id IN (' m_ids_keep_string ' )'];
% [qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% if ~isempty(emsg)
% 	disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
% else
% 	disp(['Data retrieved! In ' benrighttime(toc)]);
% end
% 
% 
% % Convert mysql NULLs to MATLAB NaNs
% empties = find(cellfun(@isempty,qrc));
% qrc(empties) = {NaN};
% Gts_id = qrc{:,1};
% Gm_id = qrc{:,2};
% Goutput = qrc{:,3};
% Gcalctime = qrc{:,4};
% Gquality = qrc{:,5};
% 
% for i=1:nts
% 	i1 = find(Gts_id==ts_ids_keep(i)); % this time series
% 	[Gm_id_S i2] = sort(Gm_id(i1)); % operations for this time series, sorted
% 	if any(Gm_id_S'-m_ids_keep) % check matches m_ids_keep
% 		disp(['Error']);
% 	end
% 	i3 = i1(i2);
% 	
% 	TS_loc(i,:) = horzcat(Goutput(i3));
% 	TS_loc_ct(i,:) = horzcat(Gcalctime(i3));
% 	TS_loc_q(i,:) = horzcat(Gquality(i3));
% 	disp([num2str(i) ' -- this took ' benrighttime(toc)]);
% end


%% METHOD 3: Get all and then filter
% % we have ts_ids_keep and m_ids_keep
% % To put in a matrix with rows (time series) and columns (metrics)
% disp('Filling local matricies TS_loc, TS_loc_ct, TS_loc_q');
% nts = length(ts_ids_keep); % number of time series
% disp([num2str(nts) ' timeseries']);
% nm = length(m_ids_keep); % number of metrics
% disp([num2str(nm) ' operations']);
% ts_ids_keep_string = bencat(ts_ids_keep,',');
% m_ids_keep_string = bencat(m_ids_keep,',');
% 
% TS_loc = zeros(nts,nm); % outputs
% TS_loc_ct = zeros(nts,nm); % calculation times
% TS_loc_q = zeros(nts,nm); % output quality label
% 
% % Retieve all data in one go from database
% disp(['Retrieving data from database...']); tic
% selectstring = ['SELECT ts_id, m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ')'];
% [qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% if ~isempty(emsg)
% 	disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
% else
% 	disp(['Data retrieved! In ' benrighttime(toc)]);
% end
% 
% % Convert mysql NULLs to MATLAB NaNs
% empties = find(cellfun(@isempty,qrc));
% qrc(empties) = {NaN};
% 
% % map to meaningful names
% Gts_id = qrc{:,1};
% Gm_id = qrc{:,2};
% Goutput = qrc{:,3};
% Gcalctime = qrc{:,4};
% Gquality = qrc{:,5};
% 
% % filter m_ids manually
% rkeep = find(ismember(m_ids_keep,Gm_id));
% Gts_id = Gts_id(rkeep);
% Gm_id = Gm_id(rkeep);
% Goutput = Goutput(rkeep);
% Gcalctime = Gcalctime(rkeep);
% Gquality = Gquality(rkeep);
% 
% 
% for i=1:nts
% 	i1 = find(Gts_id==ts_ids_keep(i)); % this time series
% 	[Gm_id_S i2] = sort(Gm_id(i1)); % operations for this time series, sorted
% 	if any(Gm_id_S'-m_ids_keep) % check matches m_ids_keep
% 		disp(['Error']);
% 	end
% 	i3 = i1(i2);
% 	
% 	TS_loc(i,:) = horzcat(Goutput(i3));
% 	TS_loc_ct(i,:) = horzcat(Gcalctime(i3));
% 	TS_loc_q(i,:) = horzcat(Gquality(i3));
% 	disp([num2str(i) ' -- this took ' benrighttime(toc)]);
% end


%% METHOD 4: row-by-row, but filter metrics in MATLAB
% dorows=1;
% if dorows
% 	times = zeros(nts,1);
% 	for i=1:nts
% 		tic
% 		selectstring = ['SELECT m_id, Output, CalculationTime, Quality FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i))];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% 		if ~isempty(emsg)
% 			disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
% 		end
% 
% 		% Convert empty entries to NaNs
% 		empties = find(cellfun(@isempty,qrc));
% 		if ~isempty(empties)
% 			qrc(empties) = {NaN};
% 		end
% 		
% 		% map to meaningful names
% 		Gm_id = qrc{:,1};
% 		Goutput = qrc{:,2};
% 		Gcalctime = qrc{:,3};
% 		Gquality = qrc{:,4};
% 		
% 		% filter m_ids manually
% 		okeep = ismember(Gm_id,m_ids_keep);
% 		Goutput = Goutput(okeep);
% 		Gcalctime = Gcalctime(okeep);
% 		Gquality = Gquality(okeep);
% 	
% 		TS_loc(i,:) = horzcat(Goutput);
% 		TS_loc_ct(i,:) = horzcat(Gcalctime);
% 		TS_loc_q(i,:) = horzcat(Gquality);
% 
% 		times(i) = toc;
% 	
% 		if mod(i,floor(nts/10))==0
% 			disp(['Approximately ' benrighttime(mean(times(1:i))*(nts-i)) ' remaining!']);
% 		end
% 	end
% 
% else % fill a column at a time
% 	times = zeros(nm,1);
% 	for i=1:nm
% 		tic
% 		selectstring = ['SELECT Output, CalculationTime, Quality FROM Results WHERE m_id = ' num2str(m_ids_keep(i)) ...
% 							' AND ts_id IN (' ts_ids_keep_string ' )'];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% 		if ~isempty(emsg)
% 			disp(['Error retrieving outputs from STORE at ' num2str(m_ids_keep(i)) ' -- Exiting']); keyboard
% 		end
% 		
% 		% Convert empty entries to NaNs
% 		empties = find(cellfun(@isempty,qrc));
% 		if ~isempty(empties)
% 			qrc(empties) = {NaN};
% 		end
% 	
% 		TS_loc(:,i) = vertcat(qrc{:,1});
% 		TS_loc_ct(:,i) = vertcat(qrc{:,2});
% 		TS_loc_q(:,i) = vertcat(qrc{:,3});
% 
% 		times(i) = toc;
% 	
% 		if mod(i,floor(nm/10))==0
% 			disp(['Approximately ' benrighttime(mean(times(1:i))*(nm-i)) ' remaining!']);
% 		end
% 	end
% 	
% 	
% end

%% Fill out guides
% First check that haven't been filtered to nothing
if nts==0 || nm==0
	disp('Enough of this, we''ve filtered ourself out of existence. That''d be right. I''m leaving.');
	SQL_closedatabase(dbc)
	return
end

% Open SQL connection
dbc = SQL_opendatabase(dbname,0);

% Retrieve Time Series Metadata
selectstring = ['SELECT FileName, Keywords, Length FROM TimeSeries WHERE ts_id IN (' ts_ids_keep_string ')'];
[tsinfo,~,~,emsg] = mysql_dbquery(dbc,selectstring);
tsf = tsinfo(:,1); % timeseries filenames
tskw = tsinfo(:,2); % timeseries keywords
tsl = vertcat(tsinfo{:,3}); % timeseries lengths

% Retrieve Operations Metadata
selectstring = ['SELECT Code, OpName, Keywords, Pointer FROM Operations WHERE m_id IN (' m_ids_keep_string ')'];
[minfo,~,~,emsg] = mysql_dbquery(dbc,selectstring);
mcode = minfo(:,1); % metric filenames
mlab = minfo(:,2); % metric labels
mkw = minfo(:,3); % metric keywords
mpoint = vertcat(minfo{:,4}); % pointer? -- [redundant if you have mlink]

% Now get master info
% (i) Which masters are implicated?
selectstring = ['SELECT Master_id, MasterLabel, MasterCode FROM MasterOperations WHERE Master_id IN ' ...
				'(SELECT DISTINCT Master_id FROM MasterPointerRelate WHERE m_id IN (' m_ids_keep_string '))'];
[masterinfo,~,~,emsg] = mysql_dbquery(dbc,selectstring);
if ~isempty(emsg)
    disp('Error retrieving Master information'); keyboard
else
    if ~isempty(masterinfo) % there are masters in out midst
        Mmid = vertcat(masterinfo{:,1});
        Mmlab = masterinfo(:,2);
        Mmcode = masterinfo(:,3);
    else
        % no error, but no master functions implicated in this subset
        Mmid = [];
        Mmlab = {};
        Mmcode = {};
    end
end

% (ii) Get master link information
selectstring = ['SELECT Master_id FROM (SELECT m_id FROM Operations WHERE m_id IN (' m_ids_keep_string ')) AS T1 ' ...
					'LEFT JOIN MasterPointerRelate ON T1.m_id = MasterPointerRelate.m_id'];
[masterlink,~,~,emsg] = mysql_dbquery(dbc,selectstring);
empties = find(cellfun(@isempty,masterlink));
if ~isempty(empties), masterlink(empties) = {0}; end
% for j=1:length(empties)
% 	masterlink{empties(j)} = 0; % means not a pointer
% end
mlink = vertcat(masterlink{:});

% (iii) renumber mlink to refer to elements of Mm rather than the indicies in Master index system
nM = length(Mmid);
rch = cell(nM,1);
for i=1:nM % number of master metrics
	rch{i} = find(mlink==Mmid(i));
end
for i=1:nM
	mlink(rch{i}) = i; % rename with index rather than the Master_id
end

% Close database
SQL_closedatabase(dbc) % close connections

%% Write output as files
% rows and columsn in TS_loc can be mapped back to STORE through tsmap and mmap
% e.g.: STORE_data(tsmap,mmap)=TS_loc, STORE_tsf(tsmap)=tsf, STORE_tsl(mmap)=mf, STORE_cts(tsmap,mmap)=TS_loc_ct

% This bit could be at the top, but actually I never use it.
if tolog
	fn = ['TS_prepared_' datestr(now,30) '.log'];
	flog = fopen(fn,'w','n');
	disp(['Logging progress to ' fn]);
	fprintf(flog,'%s\n','Welcome! This is your friendly TS_prepared log');
	fprintf(flog,'%s\n',['Subsetting and checking performed at ' datestr(now)]);
end

save('TS_loc.mat','TS_loc','-v7.3')
disp('Local version of the data: TS_loc saved');
if tolog, fprintf(flog, '%s\n', 'Local version of the data: TS_loc saved'); end

save('TS_loc_ct.mat','TS_loc_ct','-v7.3')
disp('Local version of the calc times: TS_loc_ct saved');
if tolog, fprintf(flog, '%s\n', 'Local version of TS_loc_ct saved'); end
	
save('TS_loc_q.mat','TS_loc_q','-v7.3')
disp('Local version of the output quality labels: TS_loc_q saved');
if tolog, fprintf(flog, '%s\n', 'Local version of TS_loc_q saved'); end

save('TS_loc_guides.mat','m_ids_keep','ts_ids_keep','tsf','tskw','tsl','mcode','mlab','mkw','mpoint','mlink','Mmid','Mmlab','Mmcode','-v7.3');
disp('Local version of guides saved!');
if tolog, fprintf(flog, '%s\n', 'Local version of guides saved'); end

if tolog, fclose(flog); end

%% Clean Up
% Write how many entries are to be calculated
tocalculate = sum(isnan(TS_loc_q(:)) | TS_loc_q(:)==1);
disp(['There are ' num2str(tocalculate) ' entries (= ' num2str(round(tocalculate/nm/nts*10000)/100) ' %) to calculate in TS_loc']);


% reply = input(['Shall I move on to TSQ_brawn now to calculate and then write back? [yes, ''n'' to stop]'],'s');
if ~isempty(forcalc)
	disp(['Moving on to TSQ_brawn now to calculate and then write back. No, this is not an option, this is an order.']);
	TSQ_brawn(brawninputs(1),brawninputs(2),dbname);
end


end