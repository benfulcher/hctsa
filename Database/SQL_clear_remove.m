function SQL_clear_remove(tsOrOps,idRange,doRemove,doLog)
% SQL_clear_remove
%
% Either clears results or removes entirely a given set of ts_ids
% or op_ids from the database.
%
%---INPUTS:
% tsOrOps -- either 'ts' or 'ops' for whether to eliminate either a time series
%            of a operation, respectively
% idRange -- a vector of the ts_ids or op_ids in the database to remove
% doRemove -- whether to remove entries (specify 1), or just clear their data (specify 0)
% doLog -- generate a .log file describing what was done
%
% *** Clear *** (doRemove = 0):
% The results of a particular operation or time series in the Results Table  are
% converted back to NULL. Clears *all* results from a given set of operations,
% or a given set of time series.
% This should be done whenever a time series data file is changed (or new
% results will be imcomparable to existing) or whenever a piece of code is
% altered (or its new results will be incomarable to existing results) or if
% some problem occurs.
% Note that it's better practice to clear results from all pointers to a master
% operation when the master operation is changed (especially when the master
% code is stochastic)
%
% *** Remove *** (doRemove = 1):
% Removes COMPLETELY the selected ts_ids or op_ids from the Database.

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
%% Preliminaries and input checking
% ------------------------------------------------------------------------------

if nargin < 1
	error('You must provide inputs')
end

switch tsOrOps
case 'ts'
    theWhat = 'time series';
    theid = 'ts_id';
    theTable = 'TimeSeries';
case 'ops'
    theWhat = 'operations';
    theid = 'op_id';
    theTable = 'Operations';
otherwise
    error('First input must be either ''ts'' or ''ops''')
end

% Must specify a set of time series
if nargin < 2 || min(size(idRange)) ~= 1
	fprintf(1,'Error with second input. Exiting.\n');
	return
end

if nargin < 3 % doRemove
    error('Please specify whether to remove the %s (1) or just clear their results data (0)',theWhat)
end

% Write a .log file of the clearing process by default
if nargin < 4 || isempty(doLog)
	doLog = 0;
end

% Open connection to database
[dbc, dbName] = SQL_opendatabase();

% ------------------------------------------------------------------------------
%% Provide some user feedback
% ------------------------------------------------------------------------------
if (doRemove == 0) % clear data
    input(sprintf(['Preparing to clear data for %u %s from %s.\n' ...
                                '[press any key to continue]'], ...
                                    length(idRange),theWhat,dbName),'s');
    doWhat = 'clear data';
elseif doRemove == 1
    input(sprintf(['Preparing to REMOVE %u %s from %s -- DRASTIC STUFF! ' ...
                                'I HOPE THIS IS OK?!\n[press any key to continue]'], ...
                                length(idRange),theWhat,dbName),'s');
    doWhat = 'REMOVE';
else
    error('Third input must be (0 to clear), or (1 to remove)')
end

% ------------------------------------------------------------------------------
%% Check what to clear/remove
% ------------------------------------------------------------------------------
selectString = sprintf('SELECT %s, Name FROM %s WHERE %s IN (%s)', ...
                                theid,theTable,theid,BF_cat(idRange,','));
[toDump,emsg] = mysql_dbquery(dbc,selectString);

if ~isempty(emsg)
	error('Error retrieving selected %s indices (%s) from the %s table of %s', ...
                                    	theWhat,theid,theTable,dbName)
end

if isempty(toDump)
    fprintf(1,'No %s found in the given range of %s.\n',theWhat,theid);
    return
end

toDump_id = [toDump{:,1}];
toDump_name = toDump(:,2);

input(sprintf(['About to %s %u %s stored in the Results table of %s.\n' ...
        '[press any key to show them]'],doWhat,length(toDump_name),theWhat,dbName),'s');

% ------------------------------------------------------------------------------
%% List all items to screen
% ------------------------------------------------------------------------------
for i = 1:length(toDump_id)
    fprintf(1,'[%s = %u] %s\n',theid,toDump_id(i),toDump_name{i});
end

reply = input(sprintf(['Does this look right? Check carefully -- this %s operation cannot ' ...
                            'be undone.\nType ''y'' to continue...'],doWhat),'s');
if ~strcmp(reply,'y')
	fprintf(1,'Better to be safe than sorry. Check again and come back later.\n');
	return
end

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------

if doRemove
    % ---DELETE MODE---

    % Before delete them, first get keyword information, and information about masters (for operations)
    %<><>><><><><><><>

    % First you want to get the masters
    if strcmp(tsOrOps,'ops')
        selectString = sprintf('SELECT mop_id FROM %s WHERE %s IN (%s)',theTable,theid,BF_cat(idRange,','));
        mop_ids = mysql_dbquery(dbc, selectString);
        mop_ids = mop_ids{:};
    end

	deleteString = sprintf('DELETE FROM %s WHERE %s IN (%s)',theTable,theid,BF_cat(idRange,','));
    [~,emsg] = mysql_dbexecute(dbc, deleteString);
    if isempty(emsg)
        fprintf(1,'%u %s removed from %s in %s.\n',length(toDump_id),theWhat,theTable,dbName);
    end

    SQL_FlushKeywords(tsOrOps);

    if strcmp(tsOrOps,'ops')
        % --- Update the NPointTo counters of any implicated master operations:
        fprintf(1,'Updating NPointTo counters of implicated master operations...');
        updateString = sprintf(['UPDATE MasterOperations AS m SET NPointTo = ' ...
                        '(SELECT COUNT(o.mop_id) FROM Operations AS o WHERE m.mop_id = o.mop_id) ' ...
                            'WHERE m.mop_id IN (%s)'],BF_cat(mop_ids,','));
        [~,emsg] = mysql_dbexecute(dbc, updateString);
        if ~isempty(emsg)
            error('Error counting NPointTo operations\n%s\n',emsg);
        else
            fprintf(1,' Done.\n');
        end

        % --- Delete master operations that now point to zero operations
        selectString = sprintf('SELECT mop_id, MasterCode FROM MasterOperations WHERE NPointTo = 0');
        dbOutput = mysql_dbquery(dbc, selectString);
        if isempty(dbOutput)
            fprintf(1,'All master operations still have usable operations, and remain unchanged.\n');
        else
            delete_mop_ids = dbOutput{:,1};
            delete_masterCode = dbOutput(:,2);
            deleteString = sprintf('DELETE FROM MasterOperations WHERE NPointTo = 0');
            [~,emsg] = mysql_dbexecute(dbc, deleteString);
            if ~isempty(emsg)
                error('Error deleting redundant master operations');
            else
                fprintf(1,'DELETED %u master operations that are now redundant.\n',length(delete_mop_ids));
                input('[press any key to see them]','s');
                for k = 1:length(delete_mop_ids)
                    fprintf(1,'%u/%u. [mop_id = %u]: %s\n',k,length(delete_mop_ids),delete_mop_ids(k),delete_masterCode{k});
                end
            end
        end

    end

    % Update Keyword tables (should just need to update nlink, and delete keywords that are no longer used...)
    %<><>><><><><><><>
        %<><>><><><><><><>
            %<><>><><><><><><>

    % %% Re-run keyword tables
    % if strcmp(mort, 'ts')
    %     disp(['Recalculating TimeSeriesKeywords and TsKeywordsRelate in ' dbName '. Please be patient.']);
    %     SQL_update_tskw(dbName) % updates time series keywords (will be different without the deleted time series)
    % else
    %     disp(['Recalculating OperationKeywords and OpKeywordsRelate in ' dbName '. Please be patient']);
    %     SQL_update_opkw(dbName) % updates operation keywords (will be different without the deleted operations)
    %     % disp(['Recalculating links between masters and pointers']);
    %     % SQL_linkpointermaster(dbName) % update master/pointer links
    %     % SQL_masternpointto(dbName) % counts master/pointer links for MasterOperations table
    % end
else
    %% Do the clearing
    fprintf(1,'Clearing Output, QualityCode, CalculationTime columns of the Results Table of %s...\n',dbName);
    fprintf(1,'Patience...\n');

    updateString = sprintf('UPDATE Results SET Output = NULL, QualityCode = NULL, CalculationTime = NULL WHERE %s IN (%s)',theid,BF_cat(idRange,','));
    [~,emsg] = mysql_dbexecute(dbc, updateString);

    if isempty(emsg)
    	if strcmp(tsOrOps,'ts')
    		% Get number of operations to work out how many entries were cleared
    		selectString = 'SELECT COUNT(op_id) as numOps FROM Operations';
    		numOps = mysql_dbquery(dbc,selectString); numOps = numOps{1};
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(idRange),numOps,numOps*length(idRange),dbName);
    	else
    		% Get number of time series to work out how many entries were cleared
    		selectString = 'SELECT COUNT(ts_id) as numTs FROM TimeSeries';
    		numTs = mysql_dbquery(dbc,selectString); numTs = numTs{1};
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(idRange),numTs,numTs*length(idRange),dbName);
    	end
    else
    	fprintf(1,'Error clearing results from %s... This is pretty bad news....\n%s',dbName,emsg); keyboard
    end
end

% ------------------------------------------------------------------------------
%% Close connection to the mySQL database
% ------------------------------------------------------------------------------
SQL_closedatabase(dbc) % database closed

% ------------------------------------------------------------------------------
%% Write a log file of information
% ------------------------------------------------------------------------------
if doLog
    fn = 'SQL_clear.log'; % log filename
	fprintf(1,'Writing log file to ''%s''\n',fn);

	fid = fopen(fn, 'w', 'n');
	fprintf(fid,'Document created on %s\n',datestr(now));
	if strcmp(tsOrOps,'ts')
		fprintf(fid,'Cleared outputs of %u time series\n',length(idRange));
	else
		fprintf(fid,'Cleared outputs of %u operations\n',length(idRange));
	end
	for i = 1:length(toDump_id)
        fprintf(fid,'%s\n',toDump{i});
    end
	fclose(fid);

	fprintf(1,'Logged and done and dusted!!\n');
end

fprintf(1,'Glad to have been of service to you.\n');

end
