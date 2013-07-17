function SQL_clear_remove(tsorop,vin,doremove,dbname,dolog)
% Clear the results of a particular metric or time series from the store -- i.e., convert
% results back to NULL
% Clears *all* results from a given set of operations, or a given set of time series
% This should be done whenever a time series data file is changed (or new results will be imcomparable to existing)
% or whenever a piece of code is altered (or its new results will be incomarable to existing results)
% or if some problem happens.
% Note that it's better practice to remove all pointers to a master when the master is changed (especially for
% operations that have non-deterministic outputs)

% INPUTS: tsorop -- either 'ts' or 'mets' for whether to eliminate either a time series of a metric, respectively
% 		  vin -- a vector of the ts_ids or m_ids in the database to remove
% 		  dbname -- can specify a custom database else will use default database in SQL_opendatabase
% 		  dolog -- generate a .log file describing what was done (does this by default)
% 2/12/2009 Ben Fulcher. Rehauled to use mySQL database system.
% 1/12/2012 Ben Fulcher. Minor changes.


%% Introduction, check inputs
% fprintf(1,'Welcome to ''SQL_clear''\n');

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
    theid = 'm_id';
    thetable = 'Operations';
    thename = 'OpName';
otherwise
    error('First input must be either ''ts'' or ''ops''')
end

% Must specify a set of time series
if nargin < 2 || min(size(vin)) ~= 1
	fprintf(1,'Error with second input. Exiting.\n');
	return
end

if nargin < 3 % doremove
    error('You must specify whether to remove the %s or just clear their data results',thewhat)
end

% Use default database if none specified
if nargin < 4, dbname = ''; fprintf(1,'Using default database\n'); end

% Open connection to database
[dbc, dbname] = SQL_opendatabase(dbname);

% write a .log file of the clearing process by default
if nargin < 5 || isempty(dolog)
	dolog = 0;
end

% Provide some user feedback
if doremove == 0
    reply = input(sprintf('Clearing data for %u %s from %s\n',length(vin),thewhat,dbname),'s');
    dowhating = 'clearing';
    dowhat = 'clear';
elseif doremove == 1
    reply = input(sprintf('REMOVING %u %s from %s -- SURE THIS IS OK?!\n',length(vin),thewhat,dbname),'s');
    dowhating = 'removing';
    dowhat = remove;
else
    error('Third input must be (0 == clear), or (1 == remove)')
end

% Check what to clear
SelectString = sprintf('SELECT %s FROM %s WHERE %s IN (%s)',thename,thetable,theid,BF_cat(vin,','));
[todump,~,~,emsg] = mysql_dbquery(dbc,SelectString);

if ~isempty(emsg)
	error('Error retrieving selected %s indices (%s) from the %s table of %s',thewhat,theid,thetable,dbname)
end
reply = input(sprintf(['About to clear all data from %u %s stored in the Results table of ' ...
      			dbname ' [press any key to show them]'],length(vin),thewhat),'s');

for i = 1:length(todump), fprintf(1,'%s\n',todump{i}); end
reply = input('Does this look right? Check carefully -- clearing data cannot be undone? Type ''y'' to continue...','s');
if ~strcmp(reply,'y')
	fprintf(1,'Better to be safe than sorry. Check again and come back later.\n');
	return
end

if doremove
    % Before delete them, first get keyword information, and information about masters (for operations)
    %<><>><><><><><><>
    
	DeleteString = sprintf('DELETE FROM %s WHERE %s IN (%s)',thetable,theid,BF_cat(vin,','));
    [~,emsg] = mysql_dbexecute(dbc, DeleteString);
    if isempty(emsg)
        fprintf(1,'%u %s removed from %s in %s\n',length(todump),thewhat,thetable,dbname)
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
    %     SQL_update_mkw(dbname) % updates operation keywords (will be different without the deleted operations)
    %     % disp(['Recalculating links between masters and pointers']);
    %     % SQL_linkpointermaster(dbname) % update master/pointer links
    %     % SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
    % end
else
    %% Do the clearing
    fprintf(1,'Clearing Output, QualityCode, CalculationTime columns of the Results Table of %s\n',dbname)

    UpdateString = sprintf('UPDATE Results SET Output = NULL, QualityCode = NULL, CalculationTime = NULL WHERE %s IN (%s)',theid,BF_cat(vin,','));
    [~,emsg] = mysql_dbexecute(dbc, UpdateString);

    if isempty(emsg)
    	if strcmp(tsorop,'ts')
    		% Get number of operations to work out how many entries were cleared
    		SelectString = 'SELECT COUNT(m_id) as nm FROM Operations';
    		nm = mysql_dbquery(dbc,SelectString); nm = nm{1};		
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(vin),nm,nm*length(vin),dbname);
    	else
    		% Get number of time series to work out how many entries were cleared
    		SelectString = 'SELECT COUNT(ts_id) as nts FROM TimeSeries';
    		nts = mysql_dbquery(dbc,SelectString); nts = nts{1};
    		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(vin),nts,nts*length(vin),dbname);
    	end
    else
    	fprintf(1,'Error clearing results from %s... This is pretty bad news....\n%s',dbname,emsg); keyboard
    end
end

%% Close connection to the mySQL database
SQL_closedatabase(dbc) % database closed

%% Write a log file of information
if dolog
    fn = 'SQL_clear.log'; % log filename
	fprintf(1,'Writing log file to ''%s''\n',fn)

	fid = fopen(fn, 'w', 'n');
	fprintf(fid,'Document created on %s\n',datestr(now));
	if strcmp(tsorop,'ts')
		fprintf(fid,'Cleared outputs of %u time series\n',length(vin));
	else
		fprintf(fid,'Cleared outputs of %u operations\n',length(vin));
	end
	for i = 1:length(todump), fprintf(fid,'%s\n',todump{i}); end
	fclose(fid);
    
	fprintf(1,'Logged and done and dusted!!\n');
end

fprintf(1,'Glad to of been of service to you.\n');

end