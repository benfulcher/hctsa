% fill the OperationCode table

[dbc, dbname] = SQL_opendatabase;

% 1. Scrounge through Operations looking for single operations
SelectString = 'SELECT Code FROM Operations WHERE MasterLabel IS NULL'; % not pointers


% 2. Scrounge through Masters looking for unique code

% 3. Add to OperationCode table

