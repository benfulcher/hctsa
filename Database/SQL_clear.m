function TSQ_wkshop_clear(mort,vin,dbname,dolog)
% Clear the results of a particular metric or time series from the store -- i.e., convert
% results back to NULL
% Clears *all* results from a given set of operations, or a given set of time series
% This should be done whenever a time series data file is changed (or new results will be imcomparable to existing)
% or whenever a piece of code is altered (or its new results will be incomarable to existing results)
% or if some problem happens.
% Note that it's better practice to remove all pointers to a master when the master is changed (especially for
% operations that have non-deterministic outputs)

% INPUTS: mort -- either 'ts' or 'mets' for whether to eliminate either a time series of a metric, respectively
% 		  vin -- a vector of the ts_ids or m_ids in the database to remove
% 		  dbname -- can specify a custom database else will use default database in SQL_opendatabase
% 		  dolog -- generate a .log file describing what was done (does this by default)

% 2/12/2009 Ben Fulcher. Rehauled to use mySQL database system.
% 1/12/2012 Ben Fulcher. Minor changes.


%% Introduction, check inputs
% fprintf(1,'Welcome to ''SQL_clear''\n');

% First input must be 'ts' or 'mets'
if nargin < 1 || ~ismember(mort,{'ts','mets'})
	error('First input must be either ''ts'' or ''mets''.')
end

% Must specify a set of time series
if nargin < 2 || min(size(vin))~=1
	fprintf(1,'Error with second input. Exiting.\n');
	return
end

% Use default database if none specified
if nargin < 3, dbname = ''; fprintf(1,'Using default database\n'); end

% Open connection to database
[dbc, dbname] = SQL_opendatabase(dbname);

% write a .log file of the clearing process by default
if nargin < 4 || isempty(dolog)
	dolog = 0;
end

switch mort
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
end

% Print a quick bit of user information
fprintf(1,'Clearing %u %s from %s\n',length(vin),thewhat,dbname);

% Check what to clear
SelectString = sprintf('SELECT %s FROM %s WHERE %s IN (%s)',thename,thetable,theid,bencat(vin,','));
[todump,~,~,emsg] = mysql_dbquery(dbc,SelectString);

if ~isempty(emsg)
	error(sprintf('Error retrieving selected %s indices (%s) from the %s table of %s',thewhat,theid,thetable,dbname))
end
reply = input(sprintf(['About to clear all data from %u %s stored in the Results table of ' ...
      			dbname ' [press any key to show them]'],length(vin),thewhat),'s');

for i = 1:length(todump), fprintf(1,'%s\n',todump{i}); end
reply = input('Does this look right? Check carefully -- clearing data cannot be undone? Type ''y'' to continue...','s');
if ~strcmp(reply,'y')
	fprintf(1,'Better to be safe than sorry. Check again and come back later.\n');
	return
end

   
%% Do the clearing
fprintf(1,'Clearing Output, QualityCode, CalculationTime, and LastModified columns of the Results Table of %s\n',dbname)

UpdateString = sprintf('UPDATE Results SET Output = NULL, QualityCode = NULL, CalculationTime = NULL WHERE %s IN (%s)',theid,bencat(vin,','));
[~,emsg] = mysql_dbexecute(dbc, UpdateString);

if isempty(emsg)
	if strcmp(mort,'ts')
		% get number of operations to work out how many entries were cleared
		SelectString = 'SELECT COUNT(m_id) as nm FROM Operations';
		nm = mysql_dbquery(dbc,SelectString); nm = nm{1};		
		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(vin),nm,nm*length(vin),dbname);
	else
		% get number of time series to work out how many entries were cleared
		SelectString = 'SELECT COUNT(ts_id) as nts FROM TimeSeries';
		nts = mysql_dbquery(dbc,SelectString); nts = nts{1};
		fprintf(1,'Clearing Successful! I''ve just cleared %u x %u = %u entries from %s\n',length(vin),nts,nts*length(vin),dbname);
	end
else
	fprintf(1,'Oh shit. Error clearing database... You should troubleshoot this...\n'); keyboard
end


%% Close connection to the mySQL database
SQL_closedatabase(dbc) % database closed

%% Write a log file of information
if dolog
	disp('Writing a list of time series cleared here today');
	disp('Saving information to ''TSQ_wkshop_clear.log''');

	fid = fopen('TS_wkshop_clear.log', 'w', 'n');
	fprintf(fid,'%s\n',['Document created on ' datestr(now)]);
	if strcmp(mort,'ts')
		fprintf(fid,'%s\n',['Cleared outputs of ' num2str(length(vin)) ' time series']);
	else
		fprintf(fid,'%s\n',['Cleared outputs of ' num2str(length(vin)) ' operations']);
	end
	for i = 1:length(todump), fprintf(fid,'%s\n',todump{i}); end
	fclose(fid);
	disp('Logged and done and dusted!!');
end

fprintf(1,'Glad to of been of service to you. Good day.\n');

end