% SAMPLE_RUNSCRIPT_SQL  Example script for looping over a hctsa analysis when
%                           using a linked mySQL database.

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
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
writeWhat = 'null'; % Retrieve and write back only missing (NULL) entries in the database

% Calculate across the given range of ts_ids one at a time:
tsid_range = (tsid_min:tsid_max);

% ------------------------------------------------------------------------------
%% Start calculating:
% ------------------------------------------------------------------------------
% Provide a quick message:
fprintf(1,['About to calculate across %u time series (ts_ids %u--%u) and all op_ids\n'], ...
                length(tsid_range),tsid_min,tsid_max);

% Loop across time series, one at a time:
for i = 1:length(tsid_range)
	fprintf(1,'\n\n\nWe''re looking at ts_id %u and all op_ids\n\n\n',tsid_range(i));

	% Loop over:
	% (i) Running SQL_retrieve to retrieve data from the database -> HCTSA.mat
	% (ii) Using TS_compute to calculate missing entries
	% (iii) Running SQL_store to write results back into the database

    % (i) Retrieve uncomputed entries from the database
	didWrite = SQL_retrieve(tsid_range(i),'all',writeWhat);
    if didWrite % Only calculate if SQL_retrieve found time series to retrieve:
        % (ii) Compute all the missing data in the retrieved set of
        % time series and operations:
        TS_compute(doParallelize);
        % (iii) Write the results back to the database:
        SQL_store(writeWhat);
    else
        fprintf(1,'No calculation performed at ts_id = %u\n',tsid_range(i));
    end
end
