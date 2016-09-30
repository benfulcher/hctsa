function [dbc, databaseName] = SQL_opendatabase(databaseName,beVocal,useDBToolbox)
% SQL_opendatabase 		Open a connection to a mySQL database

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Default: don't display details to the prompt:

if nargin < 2 || isempty(beVocal)
	beVocal = 0; % by default, do not display the mySQL database used to the prompt
end
if nargin < 3
    useDBToolbox = 0; % the Matlab implementation is sometimes much slower
end
% ------------------------------------------------------------------------------

[hostName, default_databaseName, username, password, customPort] = SQL_ShowConnSettings(0);

% ------------------------------------------------------------------------------
if nargin < 1 || isempty(databaseName)
	databaseName = default_databaseName; % set default database
end

% ------------------------------------------------------------------------------
if beVocal
	fprintf(1,'Using database %s\n',databaseName);
    fprintf(1,['Connecting to host ''%s'', database ''%s'' (port %u) using username' ...
            ' ''%s'' and password ''%s''...'],hostName,databaseName,customPort,username,password);
end

% ------------------------------------------------------------------------------
% Open database as dbc
% ------------------------------------------------------------------------------
dbc = mysql_dbopen(hostName,databaseName,username,password,customPort,useDBToolbox);

if isempty(dbc)
	error('Failed to load SQL database');
elseif beVocal
    fprintf(1,' Connected!\n');
end

mysql_dbexecute(dbc,['USE ' databaseName]);

end
