% ------------------------------------------------------------------------------
% sample_runscript
% ------------------------------------------------------------------------------
% 
% Sample runscript for selecting a set of operations and time series to loop
% over when running highly comparative computations.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

%% Parameters for run:
parallelize = 0; % set to 1 to parallelize computations over available CPUs using Matlab's Parellel Computing Toolbox?
DoLog = 0; % set to 1 to log results to a .log file? (usually not necessary)
tsidmin = 1; % calculate from this ts_id...
tsidmax = 100; % to this ts_id
WriteWhat = 'null'; % retrieve and write back missing (NULL) entries in the database

% Set a range of time series to calculate, as tsidr
tsidr = (tsidmin:tsidmax); % calculate across the given range of ts_ids one at a time

% Retrieve a vector of op_idds to calculate subject to additional conditions
% Here we remove operations with labels 'shit', 'tisean', 'kalafutvisscher', and 'waveletTB'
opids = SQL_getids('ops',1,{},{'shit','tisean','kalafutvisscher','waveletTB'});

% Range of op_ids retrieved at each iteration:
opidr = [min(opids), max(opids)];

%% Now start calculating
% Provide a quick message:
fprintf(1,['About to calculate across ts_ids %u--%u and op_ids %u--%u over a total of '  ...
    		 '%u iterations\n'],tsidr(1),tsidr(end),opidr(1),opidr(end),length(tsidr));

for i = 1:length(tsidr) % Loop over single time series
	fprintf(1,'\n\n\nWe''re looking at ts_id %u and %u op_ids, from %u--%u\n\n\n', ...
                                	tsidr(i),length(opids),opidr(1),opidr(2))
	
	% We loop over:
	% (i) Running TSQ_prepared to retrieve data from the database -> HCTSA_loc.mat
	% (ii) Using TSQ_brawn to calculate missing entries
	% (iii) Running TSQ_agglomerate to write results back into the database

	DidWrite = TSQ_prepared(tsidr(i),opids,WriteWhat); % Collect the null entries in the database
    if DidWrite % Only calculate if TSQ_prepared found time series to retrieve:
        TSQ_brawn(DoLog,parallelize); % computes the operations and time series retrieved
        TSQ_agglomerate(WriteWhat,DoLog); % stores the results back to the database
    else
        fprintf(1,'No time series retrieved for ts_id = %u\n',tsidr(i));
    end
end