function SQL_create_db()
	% Perform the database setup
	fprintf('Let''s set up a new database\nFirst we need to know a username and password that has create database and grant privileges\n');
	fprintf('This is probably going to be the root account\n');
	hostname = input('Hostname of MySQL server: ','s');
	admin_user = input('Username: ','s');
	admin_password = input('Password (this will appear on screen!): ','s');
	try
	    dbc = mysql_dbopen(hostname,'',admin_user,admin_password);
	    if isempty(dbc)
	        error('Please check your credentials and that the MySQL server is accessible');
	    end
	catch
	    error('Could not activate the SQL connector. Did you add it to classpath.txt? Please check the documentation');
	end

	fprintf('Connection established. Creating the database for this program\n');
	dbname = input('Enter a name for the new database: ','s');
	fprintf('And we also need a NEW username and password to create a non-admin account\n');
	local_u = input('New username: ','s');
	local_p = input('New password (this will appear on screen!): ','s');

	mysql_dbexecute(dbc,sprintf('drop database if exists %s;',dbname));
	mysql_dbexecute(dbc,sprintf('create database %s;',dbname));
	mysql_dbexecute(dbc,sprintf('grant all privileges on %s.* to ''%s''@''localhost'' identified by ''%s'';',dbname,local_u,local_p));
	mysql_dbexecute(dbc,'flush privileges;');

	SQL_closedatabase(dbc) % close the database

	fid = fopen('Database/sql_settings.conf','w');
	fprintf(fid,'%s,%s,%s,%s ',hostname,dbname,local_u,local_p);
	fclose(fid);

	try
		dbc = SQL_opendatabase;
		SQL_closedatabase(dbc);
	catch
		fprintf(1,'Error: Unable to open database after creation');
	end

