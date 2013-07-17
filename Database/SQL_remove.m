function SQL_remove(mort,vin,dbname,dolog)
% Completely removes/deletes a given set of operations or time series from the database
% The entries in the TimeSeries or Operations table is removed, and all associated entries in the
% Results table are deleted as well.

% INPUTS: mort -- either 'ts' or 'mets' for whether to eliminate either a time series of a metric, respectively
% 		  vin -- a vector of the ts_id or m_id in STORE_GUIDE/STORE_DAT to remove
% [OLD: Used to be able to use TS_sub_cntoi(tsormet,thelabs) to get the indicies of particular names of metrics/time series]

% 2/12/2009 Ben Fulcher. Rehauled to use mySQL database system.
% 1/12/2012 Ben Fulcher. A few tweaks.
error('don''t use this code')

%% Introduction
disp(['Welcome to >> TSQ_wkshop_remove << My name''s Walter Angus Alexander and I''ll be your' ...
	  ' internal monalogue for the next little while.']);
disp('I hope you enjoy your stay.');

if strcmp(mort,'ts')
	disp(['Deleting ' num2str(length(vin)) ' time series']);
elseif strcmp(mort,'mets')
	disp(['Deleting ' num2str(length(vin)) ' operations']);
else
	error('First input must be either ''ts'' or ''mets''');
end

if nargin < 2 || min(size(vin))~=1
	error('Error with second input. Exiting.');
end

if nargin < 3
	dbname = ''; % use default database by default
end

if nargin < 4 || isempty(dolog)
	dolog = 1; % logs to file by default
end

%% Open connection to MySQL Database
[dbc,dbname] = SQL_opendatabase(dbname);
disp(['Using database ' dbname]);

% Print a quick bit of user info:
if strcmp(mort,'ts')
	disp(['DELETING ' num2str(length(vin)) ' time series from ' dbname]);
elseif strcmp(mort,'mets')
	disp(['DELETING ' num2str(length(vin)) ' operations from ' dbname]);
end

%% Check what the ones are to remove
if strcmp(mort,'ts')
	SelectString = ['SELECT ts_id, Filename FROM TimeSeries WHERE ts_id IN (' BF_cat(vin,',') ')'];
	[todump,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
	input(['About to remove all trace of ' num2str(length(vin)) ' time series ' ...
				'from ' dbname ' [let''s see ''em]']);
else
	SelectString = ['SELECT m_id, OpName FROM Operations WHERE m_id IN (' BF_cat(vin,',') ')'];
	[todump,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
	input(['About to remove all trace of ' num2str(length(vin)) ' operations ' ...
				'from ' dbname ' [let''s see ''em]']);
end

theids = vertcat(todump{:,1});
todump = todump(:,2);

for i = 1:length(todump)
    disp(['[' num2str(theids(i)) ']  ' todump{i}])
end

if ~isempty(emsg)
	disp(['A problem -- maybe you''re referencing a non-existent entry?']);
end

reply = input('What do you think? Carry on? Control C now or you''re in for it...');

   
%% Do the clearing
disp('IT''S HAMMER TIME!!!');

if strcmp(mort,'ts')
	deletestring = ['DELETE FROM TimeSeries WHERE ts_id IN (' BF_cat(vin,',') ')'];
else
	deletestring = ['DELETE FROM Operations WHERE m_id IN (' BF_cat(vin,',') ')'];
end
[rs,emsg] = mysql_dbexecute(dbc, deletestring);
if ~ isempty(emsg)
	disp('Error removing entries from database... Maybe it''s a sign...? Must you always insist on destruction?'); keyboard
end
% We currently don't remove MasterOperations, we probably should check for this and update this code.
if strcmp(mort,'mets')
    fprintf('\n%s\n\n',['*** NOTE THAT WE HAVE''NT DONE ANYTHING TO THE MASTERS!!! ' ...
				'They should be evident in the MasterOperations Table from nothing pointing to them...'])
end

%% Write a log file of information
if dolog
	if strcmp(mort,'ts')
		disp('Writing a log of time series removed here today');
	else
		disp('Writing a log of operations removed here today');
	end

	fn = ['TSQ_wkshop_remove_' datestr(now,30) '.log'];
	disp(['Logging information to ''' fn '''']);

	fid = fopen(fn, 'w', 'n');
	fprintf(fid,'%s\n',['Document created on ' datestr(now)]);
	if strcmp(mort,'ts')
		fprintf(fid,'%s\n',['Cleared outputs of ' num2str(length(vin)) ' time series']);
	else
		fprintf(fid,'%s\n',['Cleared outputs of ' num2str(length(vin)) ' operations']);
	end
	for i=1:length(todump), fprintf(fid,'%s\n',todump{i}); end
	fclose(fid);
end

%% Re-run keyword tables
if strcmp(mort, 'ts')
	disp(['Recalculating TimeSeriesKeywords and TsKeywordsRelate in ' dbname '. Please be patient.']);
	SQL_update_tskw(dbname) % updates time series keywords (will be different without the deleted time series)
else
	disp(['Recalculating OperationKeywords and OpKeywordsRelate in ' dbname '. Please be patient']);
	SQL_update_mkw(dbname) % updates operation keywords (will be different without the deleted operations)
	disp(['Recalculating links between masters and pointers']);
	SQL_linkpointermaster(dbname) % update master/pointer links
	SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
end

% reply = input('Update Percentage Calculated statistics? Quite time-consuming and not really necessary... [''y'' for yes]','s');
% if strcmp(reply,'y')
% 	disp(['Doing it.']);
% 	% SQL_extra_results % recalculates the 'percentage calculated' columns in TimeSeries, TimeSeriesKeywords, Operations, and OperationKeywords
% 	SQL_fillfromresults([],[],[1 1 1],[1 1 1 1],dbname) % recalculated percentagecalculated, percentagegood, meancalctime
% end


disp('ZOMG!! DONE!!');
disp('Glad to of been of service to you. Good day.');

%% Close Database
SQL_closedatabase(dbc) % database closed

end
