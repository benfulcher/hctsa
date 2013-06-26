function SQL_create_all_tables()
% Create all the tables in the database
% Uses SQL_tablecreatestring to retrieve the appropriate mySQL CREATE TABLE statements
% Romesh Abeysuriya, March 2013
% Ben Fulcher, now uses SQL_TableCreateString, May 2013

% Specify the names of tables to create (should be valid names in SQL_TableCreateString)
TableNames = {'Operations', ...     % Operations Table
        'OperationCode', ...        % OperationCode Table
        'MasterOperations', ...     % MasterOperations Table
        'MasterPointerRelate', ...  % MasterPointerRelate Table
        'TimeSeries', ...           % TimeSeries Table
        'OperationKeywords', ...    % OperationKeywords Table
        'OpKeywordsRelate', ...     % OpKeywordsRelate Table
        'TimeSeriesKeywords', ...   % TimeSeriesKeywords
        'TsKeywordsRelate', ...     % TsKeywordsRelate Table
        'Results'};                 % Results Table

% Convert Table names to mySQL CREATE TABLE statements:
CreateString = arrayfun(@(x)SQL_TableCreateString(TableNames{x}),1:length(TableNames),'UniformOutput',0);

%% Write all of this to the database:
[dbc, dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)

fprintf(1,'Creating tables in %s\n',dbname);
nperline = 5;
for j = 1:length(CreateString)
    [rs,emsg] = mysql_dbexecute(dbc,CreateString{j});
    if ~isempty(rs)
        if j==length(CreateString)
            fprintf(1,'%s.',TableNames{j})
        else
            fprintf(1,'%s, ',TableNames{j})
        end
    else
        fprintf(1,'**** Error creating table: %s\n',TableNames{j});
        fprintf(1,'%s',emsg);
    end
    if mod(j,nperline)==0 && j < length(CreateString)
        fprintf(1,'\n'); % make new line to avoid cramping on display
    end
end
fprintf(1,'\nTables created in %s\n',dbname);

SQL_closedatabase(dbc) % close the connection to the database

end