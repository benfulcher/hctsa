% Drops and recreates the default database (set in SQL_opendatabase)
% Romesh Abeysuriya, March 2013
% Inane prompt voice added by Ben Fulcher, May 2013

[dbc, dbname] = SQL_opendatabase;
reply = input(['Are you sure you want to DELETE ALL DATA AND RESET EVERYTHING in ' dbname '?!?!?! (say ''yes'') '],'s');
if ~strcmp(reply,'yes')
    fprintf(1,'I didn''t think so... Better to be safe than sorry, hey?\n'); return
end
fprintf(1,'Omg be careful, we''re destroying everything\n');
mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s;',dbname));
fprintf(1,'%s and all the data contained within it dropped\n',dbname);
mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s;',dbname));
SQL_closedatabase(dbc) % close the database
SQL_create_all_tables;

% Add operations
SQL_add('mops','Database/INP_mops.txt','',0)
SQL_add('ops','Database/INP_ops.txt','',0)