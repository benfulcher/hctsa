fprintf(1,'This script will set up the Highly comparative time-series analysis code package from scratch\n')
fprintf(1,['In the following order, we will' ...
                '\n- Add the paths' ...
                '\n- Set up the database' ...
                '\n- Add the operations' ...
                '\n- Compile the toolboxes' ...
                '\n- Test that things are working'])

% 1. Add the paths:
fprintf(1,'Adding the paths...')
try
	startup
	fprintf('done.\n')
catch emsg
	fprintf(1,'error.\n')
    fprintf(1,'%s\n',emsg)
end

% 2. Set up the database:
fprintf(1,'Setting up the database now--NB: you need to have root access to a mySQL server to do this\n')
% Walks the user through creating the database from a root account and sets up a user account and password
SQL_create_db;
fprintf(1,['Note that if you ever want to change the database access settings, you should alter the sql_settings.conf file' ...
                ', or run SQL_create_db\n'])

% Create all tables in the database
SQL_create_all_tables;

% 3. Populate the database with operations
fprintf(1,'Populating the database with operations (please be patient)...\n')
fprintf(1,'Adding Master operations...\n'); moptic = tic;
SQL_add('mops','Database/INP_mops.txt','',0)
fprintf(1,'Master operations added in %s.\n',BF_thetime(toc(moptic)))
fprintf(1,'Adding all operations...\n'); optic = tic;
SQL_add('ops','Database/INP_ops.txt','',0)
fprintf(1,'Operations added in %s.\n',BF_thetime(toc(optic)))

% Attempt to compile the executables in Toolboxes:
fprintf(1,'Attempting to compile the binary executables needed for evaluating some operations.\n')
fprintf(1,'Please make sure that mex is set up with the right compilers for this system.\n')
fprintf(1,'Note that errors here are not the end of the world, but mean that some operations may fail to execute correctly...\n')
cd Toolboxes
compile
cd ../
fprintf(1,'Oh my goodness, everything compiled fine. The database, %s, is ready for time series to be added using SQL_add...!\n',dbname)

% Attempt to add a time series
% SQL_add('ts','INP_test_ts.txt')
