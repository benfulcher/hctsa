fprintf(1,'This script will set up the Highly Comparative Time-Series Analysis code package')
fprintf(1,['In the following order, we will ' ...
                '\n-Set up the database' ...
                '\n-Compile the toolboxes' ...
                '\n-Add the operations' ...
                '\n-Test if all is correct'])

% 1. Add the paths:
fprintf(1,'Adding the paths...')
try
	startup
	fprintf('done\n')	
catch
	fprintf(1,'Error adding the paths\n')
end


% 2. Set up the database:
fprintf(1,'Setting up the database now- you need to have root access to a MySQL server to do this\n')
SQL_create_db()
SQL_master_initiate()
cd Toolboxes
compile
cd ../

% 3. Populate database with operations
fprintf(1,'Populating the database with operations')
TSQ_add2('mets','Operations/INP_mets.txt')
