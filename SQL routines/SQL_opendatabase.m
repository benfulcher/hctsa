function [dbc dbname] = SQL_opendatabase(dbname,bevocal)
%%% SQL_opendatabase
% Opens the database as dbc for use in retrieving and storing in the mySQL
% database

if nargin<1 || isempty(dbname)
	dbname = 'tsdb_4_12_09'; % set default database
end
if nargin<2 || isempty(bevocal)
	bevocal = 1; % display the mySQL database used to the prompt
end

if bevocal
	disp(['Using database ' dbname]);
end

%% Load in database reference
dbc = mysql_dbopen('myserver.physics.ox.ac.uk',dbname,'myusername','mypassword');
if isempty(dbc)
	disp(['Error loading mySQL database']);
	return
end
mysql_dbexecute(dbc, ['use ' dbname]);

end