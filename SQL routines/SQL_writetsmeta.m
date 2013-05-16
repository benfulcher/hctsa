%%% SQL_writetsmeta
function SQL_writetsmeta(dbname)
% Writes metadata derived from reading in timeseries to the TimeSeries table in mySQL database
% Ben Fulcher 24/11/09
% Ben Fulcher 18/1/10 Added dbname argument
% Ben Fulcher 6/12/12 removed positive only field

if nargin<1
	dbname = ''; % use default
end

%% Open mySQL database
[dbc,dbname] = SQL_opendatabase(dbname); % dbc is the database

%% Do the updating
% Assumes the necessary column already exists in the database

% (*) Retrieve the filenames
disp(['Loading all entries from TimeSeries database... Please be patient...']); tic
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,'SELECT ts_id, FileName, Length from TimeSeries WHERE Length IS NULL');
nts = size(qrc,1);
disp(['Full timeseries information retrieved from ' dbname ' in ' benrighttime(toc)]);

% fill the length column of the timeseries table		
% for each time series, load the file, find the length, and write to length column of this time series
tic
for i = 1:nts
	tsid = qrc{i,1}; % unique time series ID (integer)
	tsfn = qrc{i,2}; % time series filename (string)
	tsl = qrc{i,3}; % time series length (integer/null)

	% load the time series
	x = dlmread(tsfn);

	% get information on it
	l = length(x); % length (integer)

	% update length
	if isempty(tsl) || l~=tsl % length either not yet stored, or stored incorrectly -- update
	    updatestring = ['UPDATE TimeSeries SET Length = ' num2str(l) ', LastModified = NOW() WHERE ts_id = ' num2str(tsid)];
	    [rs,emsg] = mysql_dbexecute(dbc, updatestring);
		if ~isempty(emsg)
			disp(['Error updating length for ' tsfn]);
			disp(emsg)
			SQL_closedatabase(dbc)
		end
	end
	if l~=tsl, disp(['****** Length stored incorrectly for ' tsfn ' -- fixed now']); end
	disp([num2str(i) ' -------  ' tsfn]);
end
disp(['Loaded and stored lengths in the database ' dbname ' for ' num2str(nts) ' time series!' ...
 				'\n Only took ' benrighttime(toc) ' too :)']);

%% Close database
SQL_closedatabase(dbc)