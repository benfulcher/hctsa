%% Parameters for run:
parallelize = 1; % Should I try to parallelize computations over available CPUs? (yes=1, no=0)
dolog = 0; % Should I log results to a .log file? (usually not necessary)
tslrange = [100,30000]; % set limits on the length of time series to be calculated
tsidmin = 1; % calculate from this ts_id...
tsidmax = 10; % to this ts_id

%% Open connection to default database as dbc
% [dbc,dbname] = SQL_opendatabase('');
% disp(['Reading and writing to database ' dbname]);

%% Settings for run -- how many time series / operations to retrieve at each iteration
% Set a range of time series to calculate, as tsidr
% e.g.: calculate across the first twenty ts_ids, one time series at each iteration
% tsidr = (0:1:20);
% e.g.,: calculate across the first twenty ts_ids, retrieveing two time series at each iteration
% tsidr = (0:2:20);
tsidr = ((tsidmin-1):tsidmax); % calculate across the given range of ts_ids one at a time

% % Set a range of operations to calculate, as midr
% % Retrieve all operations at each iteration:
% countstring = 'SELECT MAX(m_id) FROM Operations'; % highest m_id in Operations table
% [nm,qrf,rs,emsg] = mysql_dbquery(dbc,countstring);
% nm = nm{1}; % the highest m_id in the database
% midr = [0,nm]; % retrieves from all operations at each iteration

% retrieve a vector of m_ids to calculate subject to additional conditions
% here we remove operations with labels 'shit', 'tisean', 'kalafutvisscher', and 'waveletTB'
mids = TSQ_getids('mets',1,{},{'shit','tisean','kalafutvisscher','waveletTB'},[]); %,[midr(1)+1 midr(1+1)]);

% range of m_ids retrieved at each iteration:
midr = [min(mids),max(mids)];

% SQL_closedatabase(dbc); % database connection no longer required


%% Now start calculating
% Provide a quick message:
disp(['About to calculate across ts_ids ' num2str(tsidr(1)+1) '--' num2str(tsidr(end)) ...
		' and m_ids ' num2str(midr(1)+1) '--' num2str(midr(end)) '\n over a total of '  ...
		num2str(length(tsidr)-1) ' iterations']);

for i = 1:length(tsidr)-1 % loop over blocks of time series (tsidr)
	fprintf('\n\n\n%s\n\n\n',['We''re looking at ts_ids from ' num2str(tsidr(i)+1) ' - ' num2str(tsidr(i+1)) ...
	        ' and m_ids from ' num2str(midr(1)+1) ' - ' num2str(midr(1+1))])
	
	% retrieve a vector of ts_ids in the current range (of tsidr) with lengths between 100 and 30000
	% (but no time series labeled as 'shit' are retrieved).
	tsids = TSQ_getids('ts',tslrange,{},{'shit'},[],[tsidr(i)+1 tsidr(i+1)]);
	
	% this line uses TSQ_prepared to retrieve from the database, then runs TSQ_brawn
	% to calculate it, then runs TSQ_agglomerate to write results back to database
	TSQ_prepared(tsids,mids,0,0,'',[dolog,parallelize])
end