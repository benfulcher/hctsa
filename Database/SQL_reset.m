% Drops and recreates the default database (set in SQL_opendatabase)
% Romesh Abeysuriya, March 2013
% Inane prompt voice added by Ben Fulcher, May 2013

[dbc, dbname] = SQL_opendatabase;
reply = input(['Are you sure you want to DELETE ALL DATA AND RESET EVERYTHING in ' dbname '?!?!?! (say ''yes'')'],'s');
if ~strcmp(reply,'yes')
    disp('I didn''t think so... Better to be safe than sorry, hey?');
    return
end
disp('Omg be careful, we''re destroying everything');
mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s;',dbname));
disp([dbname ' and all the data contained within it dropped']);
mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s;',dbname));
SQL_closedatabase(dbc) % close the database
SQL_create_all_tables;
disp([dbname ' was recreated. That''s all from me ///']);