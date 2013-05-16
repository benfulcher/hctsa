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
disp(['Welcome to >> TSQ_wkshop_clear << My name''s Walter Angus Alexander and I''ll be your' ...
	  ' internal monalogue for the next little while.']);
disp('I hope you enjoy your stay.');

% First input must be 'ts' or 'mets'
if nargin<1 || ~ismember(mort,{'ts','mets'})
	error('First input must be either ''ts'' or ''mets''.')
end

% Must specify a set of time series
if nargin<2 || min(size(vin))~=1
	disp('Error with second input. Exiting.');
	return
end

% Use default database if none specified
if nargin<3, dbname = []; disp('Using default database'); end
% Open connection to database
[dbc,dbname] = SQL_opendatabase(dbname);

% write a .log file of the clearing process by default
if nargin<4 || isempty(dolog)
	dolog = 1;
end

% Print a quick bit of user information
if strcmp(mort,'ts')
	disp(['Clearing ' num2str(length(vin)) ' time series from ' dbname]);
elseif strcmp(mort,'mets')
	disp(['Clearing ' num2str(length(vin)) ' operations from ' dbname]);
end

%% Check what to clear
if strcmp(mort,'ts')
	selectstring = ['SELECT Filename FROM TimeSeries WHERE ts_id IN (' bencat(vin,',') ')'];
	[todump,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
	if ~isempty(emsg)
		SQL_closedatabase(dbc) % close connection to database
		error(['Error retrieving selected time series indices (ts_ids) from the TimeSeries table of ' dbname]);
	end
	input(['About to clear all data from ' num2str(length(vin)) ' time series stored in the Results table of ' ...
	  			dbname ' [press any key to show them]']);
else
	selectstring = ['SELECT OpName FROM Operations WHERE m_id IN (' bencat(vin,',') ')'];
	[todump,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
	if ~isempty(emsg)
		SQL_closedatabase(dbc) % close connection to database
		error(['Error retrieving selected operation indices (m_ids) from the Operations table of ' dbname]);
	end
	input(['About to clear all data from ' num2str(length(vin)) ' operations stored in the Results table of ' ...
				dbname '[press any key to show them]']);
end
char(todump)
reply = input('Does this look right? Check carefully -- clearing data cannot be undone? Type ''y'' to continue...','s');
if ~strcmp(reply,'y')
	disp('Better to be safe than sorry. Check again and come back later.');
	SQL_closedatabase(dbc) % close connection to database
	return
end

   
%% Do the clearing
disp('IT''S HAMMER TIME!!! Do not interrupt me when I''m in this mood...');
disp(['Clearing Output, Quality, CalculationTime, and LastModified columns of the Results Table of ' dbname]);

if strcmp(mort,'ts')
	updatestring = ['UPDATE Results SET Output = NULL, Quality = NULL, CalculationTime = NULL, LastModified = NOW() WHERE ts_id IN (' bencat(vin,',') ')'];
else
	updatestring = ['UPDATE Results SET Output = NULL, Quality = NULL, CalculationTime = NULL, LastModified = NOW() WHERE m_id IN (' bencat(vin,',') ')'];
end
[rs,emsg] = mysql_dbexecute(dbc, updatestring);
if isempty(emsg)
	if strcmp(mort,'ts')
		% get number of metrics to work out how many entries were cleared
		selectstring = 'SELECT COUNT(m_id) as nm FROM Operations';
		[nm,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
		nm = nm{1};		
		disp(['Clearing Successful! I''ve just cleared ' num2str(length(vin)) ' x ' num2str(nm) ' = ' num2str(nm*length(vin)) ' entries']);
	else
		% get number of time series to work out how many entries were cleared
		selectstring = 'SELECT COUNT(ts_id) as nts FROM TimeSeries';
		[nts,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
		nts = nts{1};
		disp(['Clearing Successful! I''ve just cleared ' num2str(length(vin)) ' x ' num2str(nts) ' = ' num2str(nts*length(vin)) ' entries']);
	end
else
	disp(['Oh shit. Error clearing database... You should troubleshoot this...']); keyboard
end


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

% %% Recalculate percentage calculated stats
% reply = input('Recalculate percentage calculated stats? [''y'' for yes]','s');
% if strcmp(reply,'y')
% 	SQL_fillpercentagecalc;
% end

%% Close connection to the mySQL
SQL_closedatabase(dbc) % database closed

disp('Glad to of been of service to you. Good day.');

end