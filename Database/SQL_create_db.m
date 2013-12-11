% SQL_create_db
% 
% Setup the mySQL database
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Romesh Abeysuriya
% Ben D. Fulcher <ben.d.fulcher@gmail.com>, <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function SQL_create_db()  
	fprintf(1,'Let''s set up a new database\nWe first require a username and password that has CREATE DATABASE and GRANT privileges\n');
	fprintf(1,'(This is probably going to be the root account)\n');
	hostname = input('Hostname of mySQL server (e.g., ''localhost''): ','s');
	admin_user = input('Administrator username: ','s');
	admin_password = input('Administrator password (WARNING: this will appear on screen!): ','s');
	try
	    [dbc, emsg] = mysql_dbopen(hostname,'',admin_user,admin_password);
	    if ~isempty(emsg)
	        error('Please check your username/password and that the mySQL server is accessible');
	    end
	catch
        error(['Could not activate the mySQL java connector. This must be added to Matlab''s ''javext'' directory and ' ...
                    ' the location added to Matlab''s ''classpath.txt'' file.\nCheck the documentation for details.']);
	end
    fprintf(1,'Connection to %s established\n',hostname);
    
    
    % Check it really is an account with create and grant privileges:
    % Kind of a simple one that assumes that you need "GRANT ALL PRIVILEGES"
    % Perhaps not fool-proof, so leave as a warning for now
	[a,~,~,emsg] = mysql_dbquery(dbc,'SHOW GRANTS');
    if iscell(a), a = a{1}; end % take first entry
    if ~strncmp(a,'GRANT ALL PRIVILEGES',20); % check first entry starts with "GRANT ALL PRIVILEGES"
        warning(1,'It doesn''t look like %s has administrative privileges on %s\n', admin_user, hostname);
    end
    
    % Set up a new non-admin user account for the database
	fprintf(1,'Now creating a mySQL database for highly comparative time-series analysis\n');
    fprintf(1,'Please choose a name for the database (do not include special characters).\n')
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
                'password (%s) to %s\n'],hostname,dbname,local_u,local_p,filename);
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