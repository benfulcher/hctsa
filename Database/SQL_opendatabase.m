function [dbc dbname] = SQL_opendatabase(dbname,bevocal)
%%% SQL_opendatabase
% Opens the database as dbc for use in retrieving and storing in the mySQL
% database

fid = fopen(which('sql_settings.conf'));
d = regexp(fscanf(fid,'%s'), ',', 'split');
hostname = d{1};
default_dbname = d{2};
username = d{3};
password= d{4};
fclose(fid);

if nargin<1 || isempty(dbname)
	dbname = default_dbname; % set default database
end
if nargin<2 || isempty(bevocal)
	bevocal = 0; % display the mySQL database used to the prompt
end

if bevocal
	disp(['Using database ' dbname]);
end

dbc = mysql_dbopen(hostname,dbname,username,password);

        
%% Load in database reference
if isempty(dbc)
	error(['Failed to load SQL database']);
end

mysql_dbexecute(dbc, ['use ' dbname]);

end
