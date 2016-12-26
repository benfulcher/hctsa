function out = SQL_FlushKeywords(flushWhat)
% SQL_FlushKeywords
%
% Recomputes all keywords and linkage information in the database, for either
% time series ('ts') or operations ('ops').
%
% Useful for when there's a problem with the keyword relationships (e.g., when
% an SQL_add is interrupted).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs
% ------------------------------------------------------------------------------

if nargin < 1 || isempty(flushWhat)
    flushWhat = 'ts';
end

% ------------------------------------------------------------------------------
%% Open Database
% ------------------------------------------------------------------------------
[dbc, databaseName] = SQL_opendatabase();

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

switch flushWhat
    case 'ts'
        % Check that the time series table exists in the database:
        theWhat = 'time series';
        theid = 'ts_id';
        thekid = 'tskw_id';
        theTable = 'TimeSeries';
        theKeywordTable = 'TimeSeriesKeywords';
        theRelTable = 'TsKeywordsRelate';
    case 'ops'
        theWhat = 'operations';
        theid = 'op_id';
        thekid = 'opkw_id';
        theTable = 'Operations';
        theKeywordTable = 'OperationKeywords';
        theRelTable = 'OpKeywordsRelate';
    otherwise
        error('Unknown object to update: ''%s''',flushWhat);
end


% ------------------------------------------------------------------------------
% Retrieve all keywords from the database
% ------------------------------------------------------------------------------
selectString = sprintf('SELECT %s, Keywords FROM %s',theid,theTable);
dbOutput = mysql_dbquery(dbc,selectString);

% Is the table empty?:
if isempty(dbOutput)
    % If there are entries in the keyword table, they're redundant (e.g., from an
    % SQL_clear_remove with nothing remaining afterwards)
    % You still need to drop the rel table first because of the foreign key
    % constraint (even if it's an empty table)
    warning('No entries remaining in %s -- removing all keywords.',theTable);
    theTables = mysql_dbquery(dbc,'SHOW TABLES');
    dropTable(theRelTable,theTables);
    dropTable(theKeywordTable,theTables);
    createTable(theKeywordTable);
    createTable(theRelTable);
    return
end

db_ids = [dbOutput{:,1}];
db_keywords = dbOutput(:,2);

% ------------------------------------------------------------------------------
% Split into each individual keyword
% ------------------------------------------------------------------------------
kwSplit = regexp(db_keywords,',','split','ignorecase'); % split by commas
[ukws,~,ic] = unique(horzcat(kwSplit{:})); % unique keywords
numKeywords = length(ukws); % The number of unique keywords in the new set of time series
keywordCounts = arrayfun(@(x)sum(ic==x),1:numKeywords);

fprintf(1,'I found %u unique keywords in the database of %s...\n',...
                numKeywords,theWhat);

% ------------------------------------------------------------------------------
% (if they exist) drop and then recreate the relationship and keywords tables
% (or just create them if the don't exist):
% ------------------------------------------------------------------------------
theTables = mysql_dbquery(dbc,'SHOW TABLES');
dropTable(theRelTable,theTables);
dropTable(theKeywordTable,theTables);
createTable(theKeywordTable);
createTable(theRelTable);
% (has to be in this order to avoid foreign key problems...)

% ------------------------------------------------------------------------------
% Add the unique keywords to the keywords table
% ------------------------------------------------------------------------------
insertString = sprintf('INSERT INTO %s (Keyword,NumOccur) VALUES',theKeywordTable);
toAdd = cell(numKeywords,1);
for k = 1:numKeywords
    toAdd{k} = sprintf('(''%s'',%u)',ukws{k},keywordCounts(k));
end

SQL_add_chunked(dbc,insertString,toAdd);

fprintf(1,'Added %u new keywords to %s!\n',numKeywords,theKeywordTable);

% ------------------------------------------------------------------------------
%% Fill new keyword relationships
% ------------------------------------------------------------------------------
fprintf(1,'Writing new keyword relationships to the %s table in %s...', ...
                                        theRelTable,databaseName);

relTimer = tic;

% keyword IDs in the table will be 1, ..., numKeywords since we just
% recreated the tables:
db_kwids = 1:numKeywords;

% Add for each keyword in each time series
addCell = cell(length(horzcat(kwSplit{:})),1);
k = 1;
for i = 1:length(kwSplit)
    for j = 1:length(kwSplit{i})
        addCell{k} = sprintf('(%u,%u)',db_ids(i),db_kwids(strcmp(kwSplit{i}{j},ukws)));
        k = k + 1;
    end
end

% Add them all in chunks
SQL_add_chunked(dbc,sprintf('INSERT INTO %s (%s,%s) VALUES',theRelTable,theid,thekid),addCell,500);

fprintf(1,' done in %s.\n',BF_thetime(toc(relTimer)));

% ------------------------------------------------------------------------------
function dropTable(whatTable,theTables)
    % Drops a given table (if it exists)
    if ismember(whatTable,theTables)
    	% Already exists -- drop and recreate
    	mysql_dbexecute(dbc,sprintf('DROP TABLE %s',whatTable));
    	fprintf(1,'%s table dropped.\n',whatTable);
    end
end
% ------------------------------------------------------------------------------
function createTable(whatTable)
    % Create:
    createString = SQL_TableCreateString(whatTable);
    mysql_dbexecute(dbc, createString);
    fprintf(1,'%s table created.\n',whatTable);
end
% ------------------------------------------------------------------------------

end
