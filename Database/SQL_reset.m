% ------------------------------------------------------------------------------
% SQL_reset
% ------------------------------------------------------------------------------
% Drops and recreates the default database (set in SQL_opendatabase)
% 
%---HISTORY:
% Initial idea to add by Romesh Abeysuriya, March 2013
% Inane prompt voice added by Ben Fulcher, May 2013
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
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

[dbc, dbname] = SQL_opendatabase;

% Make doubly sure the user wants to do this!!:
reply = input(['Are you sure you want to DELETE ALL DATA AND RESET EVERYTHING in ' dbname '?!?!?! (say ''yes'') '],'s');
if ~strcmp(reply,'yes')
    fprintf(1,'I didn''t think so... Better to be safe than sorry, hey?\n'); return
end
reply = input(sprintf(['Zomg be careful, we''re destroying everything.\n' ...
                        'Confirm that you want %s to be deleted? (say ''yes'') '],dbname),'s');
if ~strcmp(reply,'yes')
    fprintf(1,'I didn''t think so... Better to be safe than sorry, hey?\n'); return
end

% ------------------------------------------------------------------------------
% Drop the database:
% ------------------------------------------------------------------------------
mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s;',dbname));
fprintf(1,'%s and all the data contained within it dropped.\n',dbname);
mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s;',dbname));
SQL_closedatabase(dbc) % Close the database
SQL_create_all_tables; % Create all basic tables required by the database

% Add operations
SQL_add('mops','Database/INP_mops.txt','',0)
SQL_add('ops','Database/INP_ops.txt','',0)
