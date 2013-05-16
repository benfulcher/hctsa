fprintf(1,'This program will set up the Time Series Analysis system')
fprintf(1,'In the following order, we will\n-Set up the database\n-Compile the toolboxes\n-Add the basic metrics\n-Test if it is all correct')

fprintf(1,'Adding the paths...')
try
	startup
	fprintf('done\n')	
catch
	fprintf(1,'Error adding the paths\n')
end
fprintf(1,'Setting up the database now- you need to have root access to a MySQL server to do this\n')
SQL_create_db()
SQL_master_initiate()
cd Toolboxes
compile
cd ../

fprintf('Loading metrics into the database')
TSQ_add2('mets','Operations/INP_mets.txt')
