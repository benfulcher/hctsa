function SQL_master_initiate()
% Create all the tables in the database
% Uses SQL_tablecreatestring to retrieve the appropriate mySQL CREATE TABLE statements
% Romesh Abeysuriya, March 2013
% Ben Fulcher, now uses SQL_TableCreateString, May 2013

% Specify the names of tables to create (should be valid names in SQL_TableCreateString)
TableNames = {'Operations', ... % 1. Operations Table
        'MasterOperations', ... % 2. MasterOperations Table
        'TimeSeries', ... % 3. TimeSeries Table
        'MasterPointerRelate', ... % 4. MasterPointerRelate Table
        'OperationKeywords', ... % 5. OperationKeywords Table
        'mkwFileRelate', ... % 6. mkwFileRelate Table
        'TimeSeriesKeywords', ... % 7. TimeSeriesKeywords
        'tskwFileRelate', ...  % 8. tskwFileRelate Table
        'Results'};  % 9. Results Table

% Convert table names to mySQL CREATE TABLE statements:
CreateString = arrayfun(@(x)SQL_TableCreateString(TableNames{x}),1:length(TableNames),'UniformOutput',0);

%% Write all of this to the database:
[dbc, dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)

fprintf(1,'Creating tables in %s\n',dbname);
for j = 1:length(CreateString)
    mysql_dbexecute(dbc,CreateString{j});
    fprintf(1,'Created table: %s\n',TableNames{j});
end
fprintf(1,'Tables created in %s\n',dbname);
    
end