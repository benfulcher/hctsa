[dbc,dbname] = SQL_opendatabase;
mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s;',dbname));
mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s;',dbname));
SQL_closedatabase(dbc) % close the database
SQL_master_initiate;

