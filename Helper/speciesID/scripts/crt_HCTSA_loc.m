% ------------------------------------------------------------------------------
% crt_HCTSA_loc
% ------------------------------------------------------------------------------
%
% create a set of HCTSA files to be sent to a computer cluster so a series of
% timesereis can be calculated on different processors at the same time
% 
%% SET RANGE OF TS_IDs TO COMPUTE:
% ------------------------------------------------------------------------------
tsid_min = 121; % Calculate from this ts_id...
tsid_max = 180; % To this ts_id
nr_proc = 20;
% ------------------------------------------------------------------------------
%% Default parameters for computation:
% ------------------------------------------------------------------------------
doParallelize = 0; % Set to 1 to parallelize computations over available CPUs using Matlab's Parellel Computing Toolbox?
doLog = 0; % Set to 1 to log results to a .log file? (usually not necessary)
%writeWhat = 'null'; % Retrieve and write back only missing (NULL) entries in the database
writeWhat = 'nullerror';
% Calculate across the given range of ts_ids one at a time:
tsid_range = (tsid_min:tsid_max);
nr_tsid = tsid_max - tsid_min + 1;
ts_id_per_proc = ceil(nr_tsid/nr_proc);

% Retrieve a vector of op_idds to calculate subject to additional conditions
% Here we remove operations with labels 'shit', 'tisean', 'kalafutvisscher', and 'waveletTB'
opids = SQL_getids('ops',1,{},{'shit','kalafutvisscher','waveletTB','fa','dfa'});

%[opids,~] =  SQL_get_fast_ops_ids(4,1);
% ------------------------------------------------------------------------------
%% Start calculating:
% ------------------------------------------------------------------------------
% Provide a quick message:
fprintf(1,['About to calculate across %u time series (ts_ids %u--%u) and %u op_ids ' ...
                    '(between %u--%u) over a total of %u iterations\n'], ...
                    length(tsid_range),tsid_min,tsid_max,length(opids),min(opids),max(opids));

% Loop across time series, one at a time:
for i=1:nr_proc-1
    tsid_proc_min = (i-1)*ts_id_per_proc + 1;
    tsid_proc_max = i*ts_id_per_proc;
   
	fprintf(1,'\n\n\nWe''re looking at ts_id %u and %u op_ids, from %u--%u\n\n\n', ...
                            	tsid_range(tsid_proc_min:tsid_proc_max),length(opids),min(opids),max(opids))
	
	% Loop over:
	% (i) Running TSQ_prepared to retrieve data from the database -> HCTSA_loc.mat

	%DidWrite = TSQ_prepared(tsid_range(i),opids,writeWhat); % Collect the null entries in the database
    
    DidWrite = TSQ_prepared(tsid_range(tsid_proc_min:tsid_proc_max),opids);
    out_file_name = sprintf('HCTSA_loc_%03d.mat',i)
    movefile('HCTSA_loc.mat', out_file_name )

end

tsid_proc_min = (nr_proc-1)*ts_id_per_proc + 1;
fprintf(1,'\n\n\nWe''re looking at ts_id %u and %u op_ids, from %u--%u\n\n\n', ...
                            tsid_range(tsid_proc_min:end),length(opids),min(opids),max(opids));
DidWrite = TSQ_prepared(tsid_range(tsid_proc_min:end),opids);
out_file_name = sprintf('HCTSA_loc_%03d.mat',nr_proc)
movefile('HCTSA_loc.mat', out_file_name )



