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

function [dbc, databaseName] = SQL_opendatabase(databaseName,beVocal,useDBToolbox)

% ------------------------------------------------------------------------------
% Default: don't display details to the prompt:

if nargin < 2 || isempty(beVocal)
	beVocal = 0; % by default, do not display the mySQL database used to the prompt
end
if nargin < 3
    useDBToolbox = 0; % the Matlab implementation is sometimes much slower
end
% ------------------------------------------------------------------------------

theConfigFile = which('sql_settings.conf');
if isempty(theConfigFile)
    % no sql_settings.conf found
    error(['No sql_settings.conf file found.\n  ' ...
    'Please create sql_settings.conf file using SQL_create_db, SQL_ChangeDatabase,' ...
        ' or generate one yourself (cf. Documentation).'])
else
    fid = fopen(theConfigFile);
end

% Read the file:
s = textscan(fid,'%s%s%s%s%u','Delimiter',',','CommentStyle','%');
fclose(fid);

% Interpret the output:
hostname = s{1}{1};
default_databaseName = s{2}{1};
username = s{3}{1};
password = s{4}{1};
customPort = s{5};

% ------------------------------------------------------------------------------
if nargin < 1 || isempty(databaseName)
	databaseName = default_databaseName; % set default database
end
if isempty(customPort) || customPort==0
    customPort = 3306; % default port
end

% ------------------------------------------------------------------------------
if beVocal
	fprintf(1,'Using database %s\n',databaseName)
    fprintf(1,['Connecting to host ''%s'', database ''%s'' (port %u) using username' ...
            ' ''%s'' and password ''%s''...'],hostname,databaseName,customPort,username,password)
end

% ------------------------------------------------------------------------------
% Open database as dbc
% ------------------------------------------------------------------------------
dbc = mysql_dbopen(hostname,databaseName,username,password,customPort,useDBToolbox);

if isempty(dbc)
	error('Failed to load SQL database');
elseif beVocal
    fprintf(1,' Connected!\n');
end

mysql_dbexecute(dbc,['USE ' databaseName]);



end