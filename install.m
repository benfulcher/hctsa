fprintf(1,'This script will set up the Highly comparative time-series analysis code package')
fprintf(1,['In the following order, we will' ...
                '\n- Add the paths' ...
                '\n- Set up the database' ...
                '\n- Compile the toolboxes' ...
                '\n- Add the operations' ...
                '\n- Test if all is correct'])

% 1. Add the paths:
fprintf(1,'Adding the paths...')
try
	startup
	fprintf('done\n')
catch
	fprintf(1,'Error adding the paths\n')
end


% 2. Set up the database:
fprintf(1,'Setting up the database now--NB: you need to have root access to a mySQL server to do this\n')
% Walks the user through creating the database from a root account and sets up a user account and password
SQL_create_db;

% Create all tables in the database
SQL_create_all_tables;

% 3. Populate the database with operations
fprintf(1,'Populating the database with operations (please be patient)...\n')
fprintf(1,'Adding Master operations...'); moptic = tic;
SQL_add('mops','Database/INP_mops.txt','',0)
fprintf(1,'added in %s.',benrighttime(moptic))
fprintf(1,'Adding all operations...'); optic = tic;
SQL_add('ops','Database/INP_ops.txt','',0)
fprintf(1,'added in %s.',benrighttime(optic))

% Attempt to compile the executables in Toolboxes:
fprintf(1,'Attempting to compile the binary executables needed for evaluating some operations.\n')
fprintf(1,'Please make sure that mex is set up with the right compilers for this system.\n')
fprintf(1,'Note that errors here are not the end of the world, but mean that some operations may fail to execute correctly...\n')
cd Toolboxes
compile
cd ../