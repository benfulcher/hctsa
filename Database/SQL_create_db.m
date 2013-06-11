function SQL_create_db()
	% Setup the mysql database
	% Romesh Abeysuriya, March 2013
	% Tweaked by Ben Fulcher, May 2013
    
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
        error(['Could not activate the SQL connector. Did you add it to classpath.txt? ' ...
                    'Please check the documentation']);
	end

	fprintf('Connection established. Creating the database for this program\n');
	dbname = input('Enter a name for the new database: ','s');
	fprintf('And we also need a NEW username and password to create a non-admin account\n');
	local_u = input('New username: ','s');
	local_p = input('New password (this will appear on screen!): ','s');

	mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s',dbname));
	mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s',dbname));
	mysql_dbexecute(dbc,sprintf('GRANT ALL PRIVILEGES ON %s.* TO ''%s''@''localhost'' identified by ''%s''',dbname,local_u,local_p));
	mysql_dbexecute(dbc,'FLUSH PRIVILEGES');

	SQL_closedatabase(dbc) % close the database

    filename = 'Database/sql_settings.conf';
    disp(['Writing hostname (' hostname '), database name (' dbname '),' ...
            ' username (' local_u '), and password (' local_p ') to ' filename]);
	fid = fopen(filename,'w');
	fprintf(fid,'%s,%s,%s,%s',hostname,dbname,local_u,local_p);
	fclose(fid);

	try
		dbc = SQL_opendatabase;
		SQL_closedatabase(dbc);
        fprintf(1,['Database ' dbname ' at ' hostname ' for ' local_u ' opens and closes no problem!!']);
        fprintf(1,['We''re good to go!! :)']);
	catch
		fprintf(1,'Error: Unable to open database after creation');
	end

end