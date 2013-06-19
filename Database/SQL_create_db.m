function SQL_create_db()
	% Setup the mysql database
	% Romesh Abeysuriya, March 2013
	% Tweaked by Ben Fulcher, May-June 2013
    
	fprintf(1,'Let''s set up a new database\nFirst we need to know a username and password that has CREATE DATABASE and GRANT privileges\n');
	fprintf(1,'This is probably going to be the root account\n');
	hostname = input('Hostname of mySQL server (e.g., ''localhost''): ','s');
	admin_user = input('Username: ','s');
	admin_password = input('Password (WARNING: this will appear on screen!): ','s');
	try
	    dbc = mysql_dbopen(hostname,'',admin_user,admin_password);
	    if isempty(dbc)
	        error('Please check your credentials and that the mySQL server is accessible');
	    end
	catch
        error(['Could not activate the SQL java connector. This must be added to Matlab''s classpath.txt. ' ...
                    'Please check the documentation.']);
	end

	fprintf(1,'Connection established. Creating a mySQL database for performing highly comparative time-series analysis\n');
    fprintf(1,'Now you need to choose a name for the database. Do not include special characters.\n')
    fprintf(1,'NB: be careful, this command will replace the named database if it exists on %s.\n',hostname)
	dbname = input('Enter a name for the new database [default: hctsa]','s');
    if isempty(dbname), dbname = 'hctsa'; end
	fprintf('And we also need a NEW username and password to create a NON-ADMIN account\n');
	local_u = input('New username: ','s');
	local_p = input('New password (this will appear on screen and be stored in file!): ','s');

	mysql_dbexecute(dbc,sprintf('DROP DATABASE IF EXISTS %s',dbname));
	mysql_dbexecute(dbc,sprintf('CREATE DATABASE %s',dbname));
	mysql_dbexecute(dbc,sprintf('GRANT ALL PRIVILEGES ON %s.* TO ''%s''@''localhost'' IDENTIFIED BY ''%s''',dbname,local_u,local_p));
	mysql_dbexecute(dbc,'FLUSH PRIVILEGES');

	SQL_closedatabase(dbc) % close the database

    filename = 'Database/sql_settings.conf';
    fprintf(1,['Writing hostname (%s), database name (%s), username (%s), and ' ...
                'password (%s) to %s'],hostname,dbname,local_u,local_p,filename);
	fid = fopen(filename,'w');
	fprintf(fid,'%s,%s,%s,%s',hostname,dbname,local_u,local_p);
	fclose(fid);

	try
		dbc = SQL_opendatabase;
		SQL_closedatabase(dbc);
        fprintf(1,'Database %s at %s for %s opens and closes no problem!!\n',dbname,hostname,local_u);
        fprintf(1,'We''re good to go!! :)\n');
	catch
		fprintf(1,'Error: Unable to open database %s after creation :(\n',dbname);
	end

end