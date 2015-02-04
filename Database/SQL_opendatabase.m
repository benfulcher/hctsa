% ------------------------------------------------------------------------------
% SQL_opendatabase
% ------------------------------------------------------------------------------
% 
% Opens the database as dbc for use in retrieving and storing in the mySQL
% database
% 
%---HISTORY
% (c) 2013
% Ben D. Fulcher <ben.d.fulcher@gmail.com>, <http://www.benfulcher.com>
% 2013: Updated to use sql_settings.conf by Romesh Abeysuriya
% 2009: Max Little
% ------------------------------------------------------------------------------
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function [dbc, dbname] = SQL_opendatabase(dbname,bevocal)

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
password = d{4};
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
            ' ''%s'' and password ''%s''...'],hostname,dbname,username,password)
end

% ------------------------------------------------------------------------------
% Open database as dbc
dbc = mysql_dbopen(hostname,dbname,username,password);

if isempty(dbc)
	error('Failed to load SQL database');
elseif bevocal
    fprintf(1,' Connected!\n');
end

mysql_dbexecute(dbc,['USE ' dbname]);



end