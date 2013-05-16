%%% SQL_initiate_ts
% Sets up the TimeSeries table in the default mySQL database
% Ben Fulcher 24/11/09

%% Open database
[dbc dbname] = SQL_opendatabase; % the default database, called dbname, is opened as dbc

%% Drop table if it already exists:
% If table already exists, delete all rows in it, otherwise create it
thetables = mysql_dbquery(dbc,'Show Tables');
if ismember('timeseries',thetables)
	% already exists -- drop and recreate it
	reply = input(['Shall I really drop and reacreate the TimeSeries table from the database ' ...
					dbname '?? This will delete everything currently in the TimeSeries table. Reply with ''y'' to do this...?'],'s');
	if ~strcmp(reply,'y')
		disp('Ok, better to be safe than sorry.')
		% get out now if you don't want want to trash your existing TimeSeries table
		beep; return
	end
	mysql_dbexecute(dbc, 'DROP TABLE TimeSeries');
	disp(['TimeSeries table dropped from ' dbname]);
end

%% Create TimeSeries table in database
% Columns:
% 1) ts_id (autoincrement positive integer identifier)
% 2) FileName (string)
% 3) Keywords (string: comma-delimited keyword metadata)
% 4) Length (positive integer)
% 5) SamplingRate (float)
% 6) LastModified (datetime)
createstring = ['CREATE TABLE TimeSeries (ts_id integer not null auto_increment, Filename varchar(255), Keywords varchar(255), ' ...
				'Length integer unsigned, SamplingRate float, LastModified datetime, PRIMARY KEY (ts_id))'];
% createstring = ['CREATE TABLE TimeSeries (ts_id integer not null auto_increment, Filename varchar(255), Keywords varchar(255), ' ...
% 				'Length integer unsigned, Positive tinyint(1), Preprocess tinyint(2) unsigned, Synthetic tinyint(1), PercentageCalculated float, ' ...
% 				'PercentageGood float, MeanCalcTime float, Source varchar(500), LastModified datetime, PRIMARY KEY (ts_id))'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if isempty(emsg)
	disp(['Created table ''TimeSeries'' in ' dbname]);
else
	disp(['Error creating ''TimeSeries'' table in ' dbname]); keyboard
end

%% Close database
SQL_closedatabase(dbc)