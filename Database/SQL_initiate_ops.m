%%% SQL_initiate_ops
% Sets up the Operations and MasterOperations tables in the default mySQL server from scratch
% (drops them and the related MasterPointerRelate table if they already exist)
% Ben Fulcher, 24/11/2009

%% Open database
[dbc, dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)

%% Drop tables if they already exist
thetables = mysql_dbquery(dbc,'SHOW TABLES'); % check which tables already exist
% If Operations table already exists, drop and recreate it; otherwise just create it
if ismember('Operations',thetables)
	% already exists -- drop and recreate it
	reply = input(['Shall I really drop and reacreate the full Operations and MasterOperations tables from the' ...
						' database ' dbname '?? This will delete everything currently in these tables. Reply with ''y'' to do this...?'],'s');
	if ~strcmp(reply,'y')
		disp('Ok, better to be safe than sorry.')
		% get out now if you don't want
		SQL_closedatabase(dbc) % close the database
		beep; return
	end
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE Operations');
	if isempty(emsg)
		disp(['Operations table successfully dropped from ' dbname]);
	else
		disp(['Error dropping Operations table from ' dbname]);
		disp(emsg); beep
		keyboard
	end
end
if ismember('MasterOperations',thetables)
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE MasterOperations');
	if isempty(emsg)
		disp(['MasterOperations table dropped from ' dbname]);
	else
		disp(['Error dropping MasterOperations table from ' dbname]);
		disp(emsg); beep
		keyboard
	end
end
if ismember('MasterPointerRelate',thetables)
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE MasterPointerRelate');
	if isempty(emsg)
		disp(['MasterPointerRelate table dropped from ' dbname]);
	else
		disp(['Error dropping MasterPointerRelate table from ' dbname]);
		disp(emsg); beep
		keyboard
	end
end

%% Set up Operations Table
createstring = ['CREATE TABLE Operations (m_id integer not null auto_increment, OpName varchar(255), Pointer tinyint(1), ' ...
				'Code varchar(255), Keywords varchar(255), Stochastic tinyint(1),  LastModified datetime, ' ...
				'PRIMARY KEY (m_id))'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if isempty(emsg)
	disp(['Created ''Operations'' table in ' dbname]);
else
	disp(['Error creating ''Operations'' table in ' dbname]);
	disp(emsg); beep
	keyboard
end

%% Set up MasterOperations Table
% Create the MasterOperations table
createstring = ['CREATE TABLE MasterOperations ' ...
				'(mop_id integer not null auto_increment, MasterLabel varchar(255), ' ...
				'MasterCode varchar(255), NPointTo integer unsigned, LastModified datetime, ' ...
				'PRIMARY KEY (mop_id))'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if isempty(emsg)
	disp(['Created table ''MasterOperations'' in ' dbname]);
else
	disp(['Error creating table ''MasterOperations'' in ' dbname]);
	disp(emsg); beep
	keyboard
end


%% Close database
SQL_closedatabase(dbc) % close the database