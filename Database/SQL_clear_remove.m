% ------------------------------------------------------------------------------
% SQL_clear_remove
% ------------------------------------------------------------------------------
% 
% Either clears results or removes entirely a given set of ts_ids
% or op_ids from the database.
%
% *** Clear ***:
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
% *** Remove ***
% Removes COMPLETELY the selected ts_ids or op_ids from the Database.
% 
%
%---INPUTS:
% tsorop -- either 'ts' or 'ops' for whether to eliminate either a time series of a metric, respectively
% idRange -- a vector of the ts_ids or op_ids in the database to remove
% dbname -- can specify a custom database else will use default database in SQL_opendatabase
% doLog -- generate a .log file describing what was done (does this by default)
% 
%---HISTORY:
% 2/12/2009 Ben Fulcher. Rehauled to use mySQL database system.
% 1/12/2012 Ben Fulcher. Minor changes.
%
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Ben D. Fulcher <ben.d.fulcher@gmail.com>
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

function SQL_clear_remove(tsorop,idRange,doRemove,dbname,doLog)

% ------------------------------------------------------------------------------
%% Preliminaries and input checking
% ------------------------------------------------------------------------------

if nargin < 1
	error('You must provide inputs')
end

switch tsorop
case 'ts'
    thewhat = 'time series';
    theid = 'ts_id';
    thetable = 'TimeSeries';
    thename = 'Filename';
case 'ops'
    thewhat = 'operations';
    theid = 'op_id';
    thetable = 'Operations';
    thename = 'OpName';
otherwise
    error('First input must be either ''ts'' or ''ops''')
end

% Must specify a set of time series
if nargin < 2 || min(size(idRange)) ~= 1
	fprintf(1,'Error with second input. Exiting.\n');
	return
end

if nargin < 3 % doRemove
    error('You must specify whether to remove the %s or just clear their data results',thewhat)
end

% Use default database if none specified
if nargin < 4, dbname = ''; fprintf(1,'Using default database\n'); end

% Open connection to database
[dbc, dbname] = SQL_opendatabase(dbname);

% write a .log file of the clearing process by default
if nargin < 5 || isempty(doLog)
	doLog = 0;
end

% ------------------------------------------------------------------------------
%% Provide some user feedback
% ------------------------------------------------------------------------------
if (doRemove == 0) % clear data
    reply = input(sprintf(['Preparing to clear data for %u %s from %s. ' ...
                                '[press any key to continue]'], ...
                                    length(idRange),thewhat,dbname),'s');
    doWhat = 'clear';
elseif doRemove == 1
    reply = input(sprintf(['Preparing to REMOVE %u %s from %s -- DRASTIC STUFF! ' ...
                                'I HOPE THIS IS OK?! [press any key to continue]\n'], ...
                                length(idRange),thewhat,dbname),'s');
    doWhat = 'remove';
else
    error('Third input must be (0 to clear), or (1 to remove)')
end

% ------------------------------------------------------------------------------
%% Check what to clear/remove
% ------------------------------------------------------------------------------
selectString = sprintf('SELECT %s FROM %s WHERE %s IN (%s)', ...
                                thename,thetable,theid,BF_cat(idRange,','));
[toDump,emsg] = mysql_dbquery(dbc,selectString);

if ~isempty(emsg)
	error('Error retrieidRangeg selected %s indices (%s) from the %s table of %s', ...
                                    	thewhat,theid,thetable,dbname)
end
reply = input(sprintf(['About to clear all data from %u %s stored in the Results table of ' ...
      			dbname ' [press any key to show them]'],length(idRange),thewhat),'s');

% ------------------------------------------------------------------------------
%% List all items to screen
% ------------------------------------------------------------------------------
for i = 1:length(toDump)
    fprintf(1,'%s\n',toDump{i});
end

reply = input(['Does this look right? Check carefully -- clearing data cannot ' ...
                            'be undone? Type ''y'' to continue...'],'s');
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
    % Before delete them, first get keyword information, and information about masters (for operations)
    %<><>><><><><><><>
    
	DeleteString = sprintf('DELETE FROM %s WHERE %s IN (%s)',thetable,theid,BF_cat(idRange,','));
    [~,emsg] = mysql_dbexecute(dbc, DeleteString);
    if isempty(emsg)
        fprintf(1,'%u %s removed from %s in %s\n',length(toDump),thewhat,thetable,dbname)
    end
    
    if strcmp(tsorop,'ops')
        % What about masters??
        % 1. Get master_ids that link to deleted operations
        %<><>><><><><><><>        
        % 2. Update their NPoint counters
        %<><>><><><><><><>
        % 3. Delete those that now point to zero operations to
        %<><>><><><><><><>
    end
    
    % Update Keyword tables (should just need to update nlink, and delete keywords that are no longer used...)
    %<><>><><><><><><>
        %<><>><><><><><><>
            %<><>><><><><><><>
    
    % %% Re-run keyword tables
    % if strcmp(mort, 'ts')
    %     disp(['Recalculating TimeSeriesKeywords and TsKeywordsRelate in ' dbname '. Please be patient.']);
    %     SQL_update_tskw(dbname) % updates time series keywords (will be different without the deleted time series)
    % else
    %     disp(['Recalculating OperationKeywords and OpKeywordsRelate in ' dbname '. Please be patient']);
    %     SQL_update_opkw(dbname) % updates operation keywords (will be different without the deleted operations)
    %     % disp(['Recalculating links between masters and pointers']);
    %     % SQL_linkpointermaster(dbname) % update master/pointer links
    %     % SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
    % end
else
    %% Do the clearing
    fprintf(1,'Clearing Output, QualityCode, CalculationTime columns of the Results Table of %s...\n',dbname)
    fprintf(1,'Patience...\n');

    UpdateString = sprintf('UPDATE Results SET Output = NULL, QualityCode = NULL, CalculationTime = NULL WHERE %s IN (%s)',theid,BF_cat(idRange,','));
    [~,emsg] = mysql_dbexecute(dbc, UpdateString);

    if isempty(emsg)
    	if strcmp(tsorop,'ts')
    		% Get number of operations to work out how many entries were cleared
    		selectString = 'SELECT COUNT(op_id) as numOps FROM Operations';
    		numOps = mysql_dbquery(dbc,selectString); numOps = numOps{1};
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(idRange),numOps,numOps*length(idRange),dbname);
    	else
    		% Get number of time series to work out how many entries were cleared
    		selectString = 'SELECT COUNT(ts_id) as numTs FROM TimeSeries';
    		numTs = mysql_dbquery(dbc,selectString); numTs = numTs{1};
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(idRange),numTs,numTs*length(idRange),dbname);
    	end
    else
    	fprintf(1,'Error clearing results from %s... This is pretty bad news....\n%s',dbname,emsg); keyboard
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
	fprintf(1,'Writing log file to ''%s''\n',fn)

	fid = fopen(fn, 'w', 'n');
	fprintf(fid,'Document created on %s\n',datestr(now));
	if strcmp(tsorop,'ts')
		fprintf(fid,'Cleared outputs of %u time series\n',length(idRange));
	else
		fprintf(fid,'Cleared outputs of %u operations\n',length(idRange));
	end
	for i = 1:length(toDump), fprintf(fid,'%s\n',toDump{i}); end
	fclose(fid);
    
	fprintf(1,'Logged and done and dusted!!\n');
end

fprintf(1,'Glad to have been of service to you.\n');

end