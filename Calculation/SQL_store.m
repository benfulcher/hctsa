function SQL_store(writeWhat,dbname)
% SQL_store 	Upload data to the mySQL database.
%
% Uploads data in the HCTSA.mat file in the current directory back into the
% mySQL database. Should be done to store the result new computations done by
% TS_compute.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 1
	writeWhat = 'null'; % 'nullerror'
    % find all nulls in the database and write over them if there are values in local files
end
if ~ismember(writeWhat,{'null','error','nullerror'})
    error('Unknown specifier ''%s''',writeWhat)
end
if nargin < 3
	dbname = '';
end

% ------------------------------------------------------------------------------
%% Open a mySQL database connection
% ------------------------------------------------------------------------------
% Much faster to use java connector than the Matlab database toolbox
% (exec commands are slow; update commands are even slower)
[dbc, dbname] = SQL_opendatabase(dbname,0,0);

% ------------------------------------------------------------------------------
%% Load local files
% ------------------------------------------------------------------------------
%% Read in information from local files
fid = 1; % haha no more logging option...
loadTimer = tic;
fprintf(fid,'Loading data from HCTSA.mat...');
load('HCTSA.mat')
fprintf(fid,' Done in %s.\n',BF_thetime(toc(loadTimer)));
clear loadTimer

% ------------------------------------------------------------------------------
% Check that this file was actually retrieved from the database (such that these
% ts_ids and op_ids represent items stored in the linked mySQL database)
% ------------------------------------------------------------------------------
if ~fromDatabase
    theLoc = which('HCTSA.mat');
    error(['It looks like:\n%s\nwas generated using TS_init and does not match up ' ...
            'with ID entries a mySQL database...'],theLoc);
end

% ------------------------------------------------------------------------------
% Check that TS_CalcTime exists
% ------------------------------------------------------------------------------
% It will not exist if you've (by default) not retrieved calculation time data
% from the database, and not generated it from the computation
if ~exist('TS_CalcTime','var')
    error('No calculation time data found in HCTSA.mat');
end

% ------------------------------------------------------------------------------
%% Preliminary definitions
% ------------------------------------------------------------------------------
numTS = length(TimeSeries); % Number of time series
numOps = length(Operations); % Number of operations
ts_id_loc = [TimeSeries.ID]; % tsids in local file
op_id_loc = [Operations.ID]; % opids in local file
ts_ids_string = BF_cat(ts_id_loc,',');
op_ids_string = BF_cat(op_id_loc,',');

% ------------------------------------------------------------------------------
%% Check that nothing has been deleted since the data was loaded...
% ------------------------------------------------------------------------------
% Time series
SelectString = sprintf('SELECT COUNT(ts_id) FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string);
numTS_db = mysql_dbquery(dbc,SelectString);
numTS_db = numTS_db{1};

% Operations
SelectString = sprintf('SELECT COUNT(op_id) FROM Operations WHERE op_id IN (%s)',op_ids_string);
numOps_db = mysql_dbquery(dbc,SelectString);
numOps_db = numOps_db{1};

if (numOps == numOps_db) && (numTS == numTS_db)
    fprintf(1,'All local time series and operation ids still exist in %s. This is good.\n',dbname);
else
    if numTS_db < numTS
        fprintf(1,'There are %u time series that no longer match the database',(numTS-numTS_db));
    end

    if numOps_db < numOps
    	fprintf(1,'There are %u operations that no longer match the database',(numOps-numOps_db));
    end
    error(['It could be dangerous to write back to a changed database. ' ...
                'You should start a SQL_retrieve from scratch.'])
end

% ------------------------------------------------------------------------------
%% Find the elements that are empty in the database
%       (and hopefully full in the local file)
% ------------------------------------------------------------------------------
% Parts of calculated subsection that are empty in storage
fprintf(1,'Retrieving %s elements from the Results table in %s...',writeWhat,dbname);

switch writeWhat
case 'null'
    % collect nulls in the database
    SelectString = sprintf(['SELECT ts_id, op_id FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND QualityCode IS NULL'],ts_ids_string,op_ids_string);
case 'nullerror'
    % Collect all NULLS and previous errors
    SelectString = sprintf(['SELECT ts_id, op_id, QualityCode FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND (QualityCode IS NULL OR QualityCode = 1)'], ...
        					ts_ids_string,op_ids_string);
case 'error'
    % Collect all previous errors (assume done a SQL_retrieve using 'error' input)
    SelectString = sprintf(['SELECT ts_id, op_id FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND QualityCode = 1'], ...
        					ts_ids_string,op_ids_string);
end

retrievalTimer = tic; % Time the retrieval (should be fast)
[qrc,emsg] = mysql_dbquery(dbc,SelectString);
if ~isempty(emsg)
    fprintf(1,'\n'); error('Error selecting %s elements from %s',writeWhat,dbname);
elseif isempty(qrc)
    fprintf(1,'\nNo %s elements in this range in the database anymore! Nothing to write.\n',writeWhat);
    SQL_closedatabase(dbc); return
else
	fprintf(1,' Retrieved %u entries in %s\n',length(qrc),BF_thetime(toc(retrievalTimer)));
end
clear retrievalTimer % Stop timing


ts_id_db = vertcat(qrc{:,1}); % ts_ids (in op_id pairs) of empty database elements in this ts_id/op_id range
op_id_db = vertcat(qrc{:,2}); % op_ids (in ts_id pairs) of empty database elements in this ts_id/op_id range
numWrite = length(ts_id_db);  % Number of database elements to attempt to write back to

% ------------------------------------------------------------------------------
% Give user feedback about what database writing will occur
% ------------------------------------------------------------------------------
switch writeWhat
case 'null'
    fprintf(1,['There are %u NULL entries in Results.\nWriting calculated ' ...
                    'elements of TS_DataMat to %s...\n'],numWrite,dbname);
case 'error'
    fprintf(1,['There are %u entries in Results (all previous errors) ' ...
                    'that are being written to %s...\n'],numWrite,dbname);
    fprintf(1,['Previous results stored as errors in the database WILL NOT ' ...
                                    'be overwritten with newer errors\n']);
case 'nullerror'
    q_db = qrc(:,3); % empties (NULL) and fatal error (1)
    q_db(cellfun(@isempty,q_db)) = {0}; % turn NULLs to 0s
    q_db = vertcat(q_db{:}); % turn cells to a numeric vector
    % so now nulls in database are labeled 0, and previous errors are labeled 1
    fprintf(1,['There are %u entries in Results (either NULL or previous errors) ' ...
                    'that are being written to %s...\n'],numWrite,dbname);
    fprintf(1,['Note that previous results stored as errors in the database WILL NOT ' ...
                                'be overwritten with newer errors\n']);
    fprintf(1,'However, NULLS will be written over with any result from the local files\n');
end

localIndex = zeros(numWrite,2);
localIndex(:,1) = arrayfun(@(x)find(ts_id_loc == x,1),ts_id_db); % Indices of rows in local file for each entry in the database
localIndex(:,2) = arrayfun(@(x)find(op_id_loc == x,1),op_id_db); % Indicies of columns in local file for each entry in the database
updateMe = zeros(numWrite,1); % Label iterations that should be written to the database

% ------------------------------------------------------------------------------
% Begin writing each element to the database, one at a time, using UPDATE commands
% ------------------------------------------------------------------------------

writeBackTimer = tic; % Time how long this takes to give user feedback
numReports = 3; % the number of time remaining updates to provide the user

for i = 1:numWrite

    % Retrieve the elements
    TS_DataMat_ij = TS_DataMat(localIndex(i,1),localIndex(i,2));
    TS_Quality_ij = TS_Quality(localIndex(i,1),localIndex(i,2));
    TS_CalcTime_ij = TS_CalcTime(localIndex(i,1),localIndex(i,2));

    switch writeWhat
    case 'null'
        if isfinite(TS_DataMat_ij)
            updateMe(i) = 1; % There is a value in TS_DataMat -- write it back to the NULL entry in the database
        end
    case 'error'
        if isfinite(TS_DataMat_ij) && TS_Quality_ij~=1
            updateMe(i) = 1; % There is a now a non-error value in TS_DataMat previously returned an error (in the database)
        end
    case 'nullerror'
        if isfinite(TS_DataMat_ij) && (q_db(i)==0 || TS_Quality_ij~=1)
            updateMe(i) = 1; % there is a value in TS_DataMat -- write it to the entry in the database
        end
		% (i) Has been calculated and a value stored in TS_DataMat (isfinite()), and
		% (ii) Either the database entry is NULL or we didn't get an error (prevents writing errors over errors)
    end

    if updateMe(i)

        if isnan(TS_CalcTime_ij) % happens when there is an error in the code
            TS_CalcTime_string = 'NULL';
        else
            TS_CalcTime_string = sprintf('%f',TS_CalcTime_ij);
        end

        % I can't see any simple way around running lots of single UPDATE commands (for each entry)
    	updateString = sprintf(['UPDATE Results SET Output = %19.17g, QualityCode = %u, CalculationTime = %s ' ...
						'WHERE ts_id = %u AND op_id = %u'],TS_DataMat_ij,TS_Quality_ij, ...
						TS_CalcTime_string,ts_id_db(i),op_id_db(i));
        [~,emsg] = mysql_dbexecute(dbc, updateString);
        if ~isempty(emsg)
            SQL_closedatabase(dbc) % close the database connection before calling the error...
        	error('Error storing (ts_id,op_id) = (%u,%u) to %s??!!\n%s\n', ...
                			[TimeSeries(localIndex(i,1)).ID],[Operations(localIndex(i,2)).ID],dbname,emsg);
        end
    end

    % Give user feedback on how long is remaining:
    if (i==50) && (i < numWrite/numReports) % give an initial estimate:
        fprintf(1,['Based on the first 50 retrievals, this is taking ' ...
                'approximately %s per entry to write to the database.\n'],BF_thetime(toc(writeBackTimer)));
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(toc(writeBackTimer)/i*(numWrite-i)));
    elseif mod(i,floor(numWrite/numReports))==0 % Give numReports more time updates...
		fprintf(1,['Approximately %s remaining -- %u (/ %u possible) entries have been'  ...
			' written to %s...\n'],BF_thetime(toc(writeBackTimer)/i*(numWrite-i)),sum(updateMe),i,dbname);
	end
end

% ------------------------------------------------------------------------------
% Finished writing back to the database!
% ------------------------------------------------------------------------------

fprintf(1,['Successfully wrote %u new calculation results ' ...
            '(/ %u) to the Results table of %s in %s.\n'],...
            sum(updateMe),numWrite,dbname,BF_thetime(toc(writeBackTimer)));

if any(~updateMe) % Some entries were not written to the database
    fprintf(1,['%u entries were not written (previously-calculated errors) and remain ' ...
                            'awaiting calculation in the database.\n'],sum(~updateMe));
end

SQL_closedatabase(dbc) % Close the database connection


end
