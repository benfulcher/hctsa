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

% ------------------------------------------------------------------------------
%% SET RANGE OF TS_IDs TO COMPUTE:
% ------------------------------------------------------------------------------
tsid_min = 1; % Calculate from this ts_id...
tsid_max = 3; % To this ts_id

% ------------------------------------------------------------------------------
%% Default parameters for computation:
% ------------------------------------------------------------------------------
doParallelize = 0; % Set to 1 to parallelize computations over available CPUs using Matlab's Parellel Computing Toolbox?
doLog = 0; % Set to 1 to log results to a .log file? (usually not necessary)
writeWhat = 'null'; % Retrieve and write back only missing (NULL) entries in the database

% Calculate across the given range of ts_ids one at a time:
tsid_range = (tsid_min:tsid_max);

% Retrieve a vector of op_idds to calculate subject to additional conditions
% Here we remove operations with labels 'shit', 'tisean', 'kalafutvisscher', and 'waveletTB'
opids = SQL_getids('ops',1,{},{'shit','kalafutvisscher','waveletTB'});

% ------------------------------------------------------------------------------
%% Start calculating:
% ------------------------------------------------------------------------------
% Provide a quick message:
fprintf(1,['About to calculate across %u time series (ts_ids %u--%u) and %u op_ids ' ...
                    '(between %u--%u) over a total of %u iterations\n'], ...
                    length(tsid_range),tsid_min,tsid_max,length(opids),min(opids),max(opids));

% Loop across time series, one at a time:
for i = 1:length(tsid_range)
	fprintf(1,'\n\n\nWe''re looking at ts_id %u and %u op_ids, from %u--%u\n\n\n', ...
                            	tsid_range(i),length(opids),min(opids),max(opids))
	
	% Loop over:
	% (i) Running TSQ_prepared to retrieve data from the database -> HCTSA_loc.mat
	% (ii) Using TSQ_brawn to calculate missing entries
	% (iii) Running TSQ_agglomerate to write results back into the database

	DidWrite = TSQ_prepared(tsid_range(i),opids,writeWhat); % Collect the null entries in the database
    if DidWrite % Only calculate if TSQ_prepared found time series to retrieve:
        TSQ_brawn(doLog,doParallelize); % Compute the operations and time series retrieved
        TSQ_agglomerate(writeWhat,doLog); % Store the results back to the database
    else
        fprintf(1,'No calculation performed at ts_id = %u\n',tsid_range(i));
    end
end