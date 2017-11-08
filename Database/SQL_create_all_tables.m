function SQL_create_all_tables()
% SQL_create_all_tables      Create all the tables in the database
%
% Uses SQL_TableCreateString to retrieve the appropriate mySQL CREATE TABLE
% statements.

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Specify the names of tables to create (should be valid names in SQL_TableCreateString)
tableNames = {'MasterOperations', ...     % MasterOperations Table
              'Operations', ...           % Operations Table
              'TimeSeries', ...           % TimeSeries Table
              'OperationKeywords', ...    % OperationKeywords Table
              'OpKeywordsRelate', ...     % OpKeywordsRelate Table
              'TimeSeriesKeywords', ...   % TimeSeriesKeywords
              'TsKeywordsRelate', ...     % TsKeywordsRelate Table
              'Results',...               % Results Table
              'GitInfo'};                 % Git repository information

% Convert Table names to mySQL CREATE TABLE statements:
createString = arrayfun(@(x)SQL_TableCreateString(tableNames{x}),...
                        1:length(tableNames),'UniformOutput',0);
existString = arrayfun(@(x)['SHOW TABLES LIKE ''' tableNames{x} ''''],...
                        1:length(tableNames),'UniformOutput',0);

% ------------------------------------------------------------------------------
%% Write all of this to the database:
% ------------------------------------------------------------------------------
[dbc,dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)

numPerLine = 3; % Make a new line after adding this many tables to the database
fprintf(1,'Creating tables in %s:\n',dbname);
for j = 1:length(createString)

    % First check whether the table already exists:
    output = mysql_dbquery(dbc,existString{j});
    if ~isempty(output) % Table already exists
        if j == length(createString)
            fprintf(1,'(%s already exists).',tableNames{j});
        else
            fprintf(1,'(%s already exists), ',tableNames{j});
        end
    else
        % Table does not yet exist; attempt to create it:
        [rs,emsg] = mysql_dbexecute(dbc,createString{j});
        if ~isempty(rs)
            if j == length(createString)
                fprintf(1,'%s.',tableNames{j});
            else
                fprintf(1,'%s, ',tableNames{j});
            end
        else
            fprintf(1,'**** Error creating table: %s [%s]\n',tableNames{j},emsg);
        end
    end

    % Make new line to avoid cramping on display
    if (mod(j,numPerLine) == 0) && (j < length(createString))
        fprintf(1,'\n');
    end
end
fprintf(1,'\nTables created in %s.\n',dbname);

%-------------------------------------------------------------------------------
% Add git repository information
%-------------------------------------------------------------------------------
gitInfo = TS_AddGitInfo();
if isempty(gitInfo)
    warning('No git information...?')
end
values = {gitInfo.branch,gitInfo.hash,gitInfo.remote,gitInfo.url};
values = cellfun(@makeSQLFriendly,values,'UniformOutput',0);
insertString = sprintf(['INSERT INTO GitInfo (branch,hash,remote,url) VALUES ',...
                        '(%s,%s,%s,%s)'],values{1},values{2},values{3},values{4});
[rs,emsg] = mysql_dbexecute(dbc,insertString);

%-------------------------------------------------------------------------------
% Close the connection to the database
SQL_closedatabase(dbc);

%-------------------------------------------------------------------------------
function sFriendly = makeSQLFriendly(s)
    % Turn an element of gitInfo into an insert-friendly thingo
    if isempty(s)
        sFriendly = 'NULL';
    elseif ischar(s)
        sFriendly = sprintf('''%s''',s);
    end
end

end
