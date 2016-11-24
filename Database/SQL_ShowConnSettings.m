function [hostName, databaseName, username, password, customPort] = SQL_ShowConnSettings(doDisplay)
% SQL_ShowConnSettings  Display the database connection settings,
% read in from the sql_settings.conf file.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs:
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(doDisplay)
    doDisplay = 1;
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
hostName = s{1}{1};
databaseName = s{2}{1};
username = s{3}{1};
password = s{4}{1};
customPort = s{5};

if isempty(customPort) || customPort==0
    customPort = 3306; % default port
end

% ------------------------------------------------------------------------------

if doDisplay
    fprintf(1,'Connecting to %s at %s as %s (through port %u)\n',...
                    databaseName,hostName,username,customPort);
end


end
