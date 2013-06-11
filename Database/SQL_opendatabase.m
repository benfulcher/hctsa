function [dbc, dbname] = SQL_opendatabase(dbname,bevocal)
% Opens the database as dbc for use in retrieving and storing in the mySQL database
% Ben Fulcher, 2009, adapting code provided by Max Little
% Updated to use sql_settings.conf by Romesh Abeysuriya

fid = fopen(which('sql_settings.conf'));
d = regexp(fscanf(fid,'%s'), ',', 'split');
hostname = d{1};
default_dbname = d{2};
username = d{3};
password= d{4};
fclose(fid);

if nargin < 1 || isempty(dbname)
	dbname = default_dbname; % set default database
end
if nargin < 2 || isempty(bevocal)
	bevocal = 0; % display the mySQL database used to the prompt
end

if bevocal
	disp(['Using database ' dbname]);
end

%% Open database as dbc
dbc = mysql_dbopen(hostname,dbname,username,password);

if isempty(dbc)
	error('Failed to load SQL database');
end

mysql_dbexecute(dbc, ['USE ' dbname]);

end