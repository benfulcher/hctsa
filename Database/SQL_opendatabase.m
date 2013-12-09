function [dbc, dbname] = SQL_opendatabase(dbname,bevocal)
% Opens the database as dbc for use in retrieving and storing in the mySQL database
% Ben Fulcher, 2009, adapting code provided by Max Little
% Updated to use sql_settings.conf by Romesh Abeysuriya

theconfigfile = which('sql_settings.conf');
if isempty(theconfigfile)
    % no sql_settings.conf found
    error('No sql_settings.conf file found. Please create sql_settings.conf file using SQL_create_db.')
else
    fid = fopen(theconfigfile);
end
s = textscan(fid,'%s%s%s%s','Delimiter',',','CommentStyle','%','CollectOutput',1); % 'HeaderLines',1,
d = s{1};
% d = regexp(s,',','split');
% d = regexp(fscanf(fid,'%s'), ',', 'split');
hostname = d{1};
default_dbname = d{2};
username = d{3};
password= d{4};
fclose(fid);

if nargin < 1 || isempty(dbname)
	dbname = default_dbname; % set default database
end
if nargin < 2 || isempty(bevocal)
	bevocal = 0; % by default, do not display the mySQL database used to the prompt
end

if bevocal
	fprintf(1,'Using database %s\n',dbname)
    fprintf(1,['Connecting to host ''%s'', database ''%s'', using username' ...
            ' ''%s'' and password ''%s''\n'],hostname,dbname,username,password)
end

%% Open database as dbc
dbc = mysql_dbopen(hostname,dbname,username,password);

if isempty(dbc)
	error('Failed to load SQL database');
end

mysql_dbexecute(dbc,['USE ' dbname]);

end