%% Parameters for run:
parallelize = 0; % set to 1 to parallelize computations over available CPUs using Matlab's Parellel Computing Toolbox?
dolog = 0; % set to 1 to log results to a .log file? (usually not necessary)
tslrange = [100, 30000]; % set limits on the length of time series to be calculated
tsidmin = 1; % calculate from this ts_id...
tsidmax = 5; % to this ts_id
writewhat = 'null'; % retrieve and write back NULL entries in the database

%% Settings for run -- how many time series / operations to retrieve at each iteration
% Set a range of time series to calculate, as tsidr
% e.g.: calculate across the first twenty ts_ids, one time series at each iteration
% tsidr = (0:1:20);
% e.g.,: calculate across the first twenty ts_ids, retrieveing two time series at each iteration
% tsidr = (0:2:20);
tsidr = ((tsidmin-1):tsidmax); % calculate across the given range of ts_ids one at a time

% % Set a range of operations to calculate, as midr

% retrieve a vector of m_ids to calculate subject to additional conditions
% here we remove operations with labels 'shit', 'tisean', 'kalafutvisscher', and 'waveletTB'
mids = TSQ_getids('mets',1,{},{'shit','tisean','kalafutvisscher','waveletTB','locdep','spreaddep'},[]);

% range of m_ids retrieved at each iteration:
midr = [min(mids), max(mids)];

%% Now start calculating
% Provide a quick message:
fprintf(1,['About to calculate across ts_ids %u--%u and m_ids %u--%u over a total of '  ...
    		 '%u iterations'],tsidr(1)+1,tsidr(end),midr(1)+1,midr(end),length(tsidr)-1);

for i = 1:length(tsidr)-1 % loop over blocks of time series (tsidr)
	fprintf(1,'\n\n\nWe''re looking at ts_ids from %u--%u and m_ids from %u--%u\n\n\n', ...
                            	tsidr(i)+1,tsidr(i+1),midr(1),midr(2))
	
	% retrieve a vector of ts_ids in the current range (of tsidr) with lengths between 100 and 30000
	% (but no time series labeled as 'shit' are retrieved).
	tsids = TSQ_getids('ts',tslrange,{},{'shit'},[],[tsidr(i)+1 tsidr(i+1)]);
	
	% this line uses TSQ_prepared to retrieve from the database, then runs TSQ_brawn
	% to calculate it, then runs TSQ_agglomerate to write results back to database

	TSQ_prepared(tsids,mids,writewhat); % Collect the null entries in the database
    TSQ_brawn; % computes the operations and time series retrieved
    TSQ_agglomerate(writewhat,dolog); % stores the results back to the database
end