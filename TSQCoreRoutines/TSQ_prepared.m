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
% 6) tolog [opt]: whether or not to log to a local file a summary of the input/output. (default = 0 -- don't log)
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
if nargin < 3 || isempty(tolog)
	tolog = 0;  % no log created
end
if nargin < 4
	forcalc = []; % retrieve full sets of things, not just the empty entries in the database
end
if nargin < 5
	dbname = []; % Use default database
end
if nargin < 6 || isempty(brawninputs)
	brawninputs = [1, 1]; % log and parallelize -- only relevant if forcalc isn't empty
end
if nargin < 7 || isempty(doinone)
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
TS_loc = zeros(nts,nm); % outputs
TS_loc_ct = zeros(nts,nm); % calculation times
TS_loc_q = zeros(nts,nm); % output quality label

% Open database connection
[dbc,dbname] = SQL_opendatabase(dbname);

% Tell me about it
fprintf(1,'We have %u time series and %u operations to retrieve from %s\n',nts,nm,dbname);
fprintf(1,'Filling and saving to local matricies TS_loc, TS_loc_ct, TS_loc_q from Results table in %s\n',dbname);

if ~isempty(forcalc)
	fprintf(1,'Retrieving elements to calculate from the database. Please be patient... I''m timing...');
	tic
	% Just retrieve bits that need to be calculated
	% May be faster to break it up and do each ts_id seperately... (this way the amount of data in each transfer is limited)
	SelectString = ['SELECT ts_id, m_id, Output, CalculationTime, QualityCode FROM Results WHERE ts_id IN (' ts_ids_keep_string ')' ...
						' AND m_id IN (' m_ids_keep_string ' ) AND QualityCode IS NULL'];
	if forcalc > 0
		SelectString = [SelectString, ' LIMIT ' num2str(forcalc)];
	end
    
	[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
	if ~isempty(emsg)
		fprintf(1,'Error retrieving outputs from database...???\n');
        fprintf(1,'%s\n',emsg)
        keyboard
	else
		ngot = size(qrc,1);
		if ngot == 0
			fprintf(1,'Nothing to retrieve?!!\n');
			SQL_closedatabase(dbc)
			return
		else
			fprintf(1,'Retrieved %u elements from %s!!! Took %s\n',ngot,dbname,benrighttime(toc));
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
	for i = 1:nts
		Gindex_ts(Gts_ids == ts_ids_keep(i)) = i; % switch from the ts_id to the index
	end
	Gindex_m = Gm_ids;
	for i = 1:nm
		Gindex_m(Gm_ids == m_ids_keep(i)) = i; % switch from the m_id to the index
	end
	% Gindex_pairs = [Gindex_ts,Gindex_m]
	
	% We have the index pairs, we now must just write these into the matricies
	% First, initiate matricies with Infs (means bogus -- unattained)
	
	fprintf(1,'Filling local files TS_loc, TS_loc_ct, and TS_loc_q...');
	tic
	TS_loc(:) = Inf; TS_loc_ct(:) = Inf; TS_loc_q = Inf;
	for i = 1:ngot
		TS_loc(Gindex_ts(i),Gindex_m(i)) = Goutput(i);
		TS_loc_ct(Gindex_ts(i),Gindex_m(i)) = Gcalctime(i);
		TS_loc_q(Gindex_ts(i),Gindex_m(i)) = Gquality(i);
	end
	fprintf(1,'filled. Took %s\n',benrighttime(toc));
	
	
	% ()()() Filter down TS_loc_ and ts_ids_keep, m_ids_keep ()()()
	fprintf(1,'Filtering unused time series and operations\n');
	tic
	% Time Series
	Utsids = unique(Gindex_ts); % rows with things to calculate in them (will be rows of Infs in TS_loc_)
	nUtsids = length(Utsids);
	if nUtsids < nts
		fprintf(1,'Cutting down from %u to %u time series with calculation potential\n',nts,nUtsids);
		ts_ids_keep = ts_ids_keep(Utsids);
		TS_loc = TS_loc(Utsids,:);
		TS_loc_ct = TS_loc_ct(Utsids,:);
		TS_loc_q = TS_loc_q(Utsids,:);
	end
	
	% Operations
	Umids = unique(Gindex_m); % rows with things to calculate in them (will be rows of Infs in TS_loc_)
	nUmids = length(Umids);
	if nUmids < nm
		fprintf(1,'Cutting down from %u to %u operations with calculation potential\n',nm,nUmids);
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
	
	fprintf(1,'Filtering complete in %s\n',benrighttime(toc));	
	
else
    % retrieve everything in the range given (not just empty bits)
	% choose whether to do by rows or columns	
	
	%% Strategy: fill the matrices one row at a time
	% Opening and closing individual connections helps to avoid a massive java heap space error... Hopefully...
	times = zeros(nts,1);
    % if doinone
        % dbc = SQL_opendatabase(dbname);
    % end
	for i = 1:nts
		tic
		SelectString = sprintf(['SELECT m_id, Output, CalculationTime, QualityCode FROM Results WHERE ' ...
                            		'ts_id = %u AND m_id IN (%s)'],ts_ids_keep(i),m_ids_keep_string);
		[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);

		if ~isempty(emsg)
			fprintf(1,'Error retrieving outputs from the Results table in %s at ts_id = %u...\n',dbname,ts_ids_keep(i));
            keyboard
		end
	
        if isempty(qrc)
            error(['No entries found in the Results table in ' dbname ...
                    '... Have you specified the right set of time series and operations???'])
        end
        
		% Convert empty entries to NaNs
		qrc(cellfun(@isempty,qrc)) = {NaN};
	
		% give sensible names
		Gm_id = horzcat(qrc{:,1});
		Goutput = horzcat(qrc{:,2});
		Gcalctime = horzcat(qrc{:,3});
		Gquality = horzcat(qrc{:,4});
	
		% make sure m_ids are in ascending order to match
		% (should be, if table is ordered properly...)
		[S_m_ids, ix] = sort(Gm_id);
	
		% Check matches
		if ~all(S_m_ids' - m_ids_keep == 0)
			fprintf(1,'Problem with %u\n',ts_ids_keep(i));
			keyboard
		end
	
		Goutput = Goutput(ix);
		Gcalctime = Gcalctime(ix);
		Gquality = Gquality(ix);

		TS_loc(i,:) = Goutput;
		TS_loc_ct(i,:) = Gcalctime;
		TS_loc_q(i,:) = Gquality;

		times(i) = toc;
		if mod(i,floor(nts/10)) == 0
			fprintf(1,'Approximately %s remaining!\n',benrighttime(mean(times(1:i))*(nts-i)));
		end
	end
    % if doinone
        % SQL_closedatabase(dbc);
    % end

	fprintf(1,'Done. Took %s altogether.\n',benrighttime(sum(times)));

	% else % fill a column at a time
	% 	times = zeros(nm,1);
	% 	for i=1:nm
	% 		tic
	% 		SelectString = ['SELECT Output, CalculationTime, QualityCode FROM Results WHERE ts_id IN (' ts_ids_keep_string ' ) ' ...
	% 						'AND m_id = ' num2str(m_ids_keep(i))];
	% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% 		SelectString = ['SELECT m_id, Output, CalculationTime, QualityCode FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i)) ...
% 							' AND m_id IN (' m_ids_keep_string ' ) AND Output IS NULL'];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% 		SelectString = ['SELECT Output, CalculationTime, QualityCode FROM Results WHERE ts_id IN (' ts_ids_keep_string ' ) ' ...
% 						'AND m_id = ' num2str(m_ids_keep(i))];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% SelectString = ['SELECT ts_id, m_id, Output, CalculationTime, QualityCode FROM Results WHERE ts_id IN (' ts_ids_keep_string ') ' ...
% 					'AND m_id IN (' m_ids_keep_string ' )'];
% [qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% SelectString = ['SELECT ts_id, m_id, Output, CalculationTime, QualityCode FROM Results WHERE ts_id IN (' ts_ids_keep_string ')'];
% [qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% 		SelectString = ['SELECT m_id, Output, CalculationTime, QualityCode FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i))];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
% 		SelectString = ['SELECT Output, CalculationTime, QualityCode FROM Results WHERE m_id = ' num2str(m_ids_keep(i)) ...
% 							' AND ts_id IN (' ts_ids_keep_string ' )'];
% 		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
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
if nts == 0 || nm == 0
	fprintf(1,'After filtering, there is nothing left to do! Enough of this then...');
	SQL_closedatabase(dbc)
	return
end

% Open a database connection
% dbc = SQL_opendatabase(dbname,0);

% Retrieve Time Series Metadata
SelectString = sprintf('SELECT FileName, Keywords, Length FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_keep_string);
[tsinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
tsf = tsinfo(:,1); % timeseries filenames
tskw = tsinfo(:,2); % timeseries keywords
tsl = vertcat(tsinfo{:,3}); % timeseries lengths

% Retrieve Operations Metadata
SelectString = sprintf('SELECT Code, OpName, Keywords FROM Operations WHERE m_id IN (%s)',m_ids_keep_string);
[minfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
mcode = minfo(:,1); % metric filenames
mlab = minfo(:,2); % metric labels
mkw = minfo(:,3); % metric keywords
% mpoint = vertcat(minfo{:,4}); % pointer? -- [redundant if you have mlink]

% Now get master info
% (i) Which masters are implicated?
SelectString = ['SELECT mop_id, MasterLabel, MasterCode FROM MasterOperations WHERE mop_id IN ' ...
				'(SELECT DISTINCT mop_id FROM MasterPointerRelate WHERE m_id IN (' m_ids_keep_string '))'];
[masterinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
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
	rch{i} = find(mlink==Mmid(i));
    mlink(rch{i}) = i; % rename with index rather than the mop_id
end

% Close database connection
SQL_closedatabase(dbc)

%% Write output as files
% rows and columsn in TS_loc can be mapped back to STORE through tsmap and mmap
% e.g.: STORE_data(tsmap,mmap)=TS_loc, STORE_tsf(tsmap)=tsf, STORE_tsl(mmap)=mf, STORE_cts(tsmap,mmap)=TS_loc_ct

% This bit could be at the top, but actually I never use it.
if tolog
	fn = sprintf('TS_prepared_%s.log',datestr(now,30)); % filename
	flog = fopen(fn,'w','n');
	fprintf(1,'Logging progress to %s\n',fn);
	fprintf(flog,'Welcome! This is your friendly TSQ_prepared log.\n');
	fprintf(flog,'Subsetting and checking performed at %s\n',datestr(now));
end

fprintf(1,'Saving local versions of the data...:\n');
save('TS_loc.mat','TS_loc','-v7.3')
fprintf(1,'TS_loc (data)')
if tolog, fprintf(flog, '%s\n', 'Local version of the data: TS_loc saved'); end

save('TS_loc_ct.mat','TS_loc_ct','-v7.3')
fprintf(1,', TS_loc_ct (calculation times)')
if tolog, fprintf(flog, '%s\n', 'Local version of TS_loc_ct saved'); end
	
save('TS_loc_q.mat','TS_loc_q','-v7.3')
fprintf(1,', TS_loc_q (quality codes)')
if tolog, fprintf(flog, '%s\n', 'Local version of TS_loc_q saved'); end

save('TS_loc_guides.mat','m_ids_keep','ts_ids_keep','tsf','tskw','tsl','mcode','mlab','mkw','mlink','Mmid','Mmlab','Mmcode','-v7.3');
fprintf(1,', TS_loc_guides (information guides).\nAll saved.\n')
if tolog, fprintf(flog, '%s\n', 'Local version of guides saved'); end

if tolog, fclose(flog); end

%% Clean Up
% Write how many entries are to be calculated
tocalculate = sum(isnan(TS_loc_q(:)) | TS_loc_q(:)==1);
fprintf(1,'There are %u entries (= %5.2f%%) to calculate in TS_loc (%ux%u)\n',tocalculate,tocalculate/nm/nts*100,nts,nm);

% reply = input(['Shall I move on to TSQ_brawn now to calculate and then write back? [yes, ''n'' to stop]'],'s');
if ~isempty(forcalc)
	fprintf(1,'Automatically moving on to TSQ_brawn now to calculate the missing entries and then write them back to the database\n');
	TSQ_brawn(brawninputs(1),brawninputs(2),dbname);
end


end