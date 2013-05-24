function SQL_update_mkw(dbname)
% Recreates the keywords tables and their links to metrics
% To be run when operations are either added or removed
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 turned into a function; added dbname input

if nargin < 1
	dbname = ''; % opens default specified in SQL_opendatabase
end

%% Open database
dbc = SQL_opendatabase(dbname); % dbc is the database

%% Drop existing tables
% (1) mkwFileRelate
[thetables,~,~,emsg] = mysql_dbquery(dbc,'SHOW TABLES');
if ismember('mkwFileRelate',thetables)
	% alread exists -- drop and recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE mkwFileRelate');
	if ~isempty(rs)
		disp('mkwFileRelate table dropped');
	else
        disp(emsg);
		error('Error dropping table mkwFileRelate');
	end
else
    disp(['No mkwFileRelate table in ' dbname]);
end

% (2) Operation Keywords
if ismember('OperationKeywords',thetables)
	% alread exists -- drop then recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE OperationKeywords');
	if ~isempty(rs)
		disp('Successfully dropped OperationKeywords table');
	else
        disp(emsg);
		error('Error dropping OperationKeywords table');
	end
else
    disp(['No OperationKeywords table in ' dbname]);
end


%% Find all unique keywords strings, split into a table of unique, separate keywords
disp('Looking for unique keywords -- to place in new table OperationKeywords');
 
selectstring = 'SELECT DISTINCT Keywords FROM Operations';
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);

splitkws = {};
k = 1;
for i = 1:length(qrc)
    kws = regexp(qrc{i},',','split','ignorecase');
    for j = 1:length(kws)
        splitkws{k} = kws{j};
        k = k + 1;
    end
end
ukws = unique(splitkws); % cell of unique keyword strings


% Create OperationKeywords Table
CreateString = SQL_TableCreateString('OperationKeywords');
[rs,emsg] = mysql_dbexecute(dbc, CreateString);
if ~isempty(rs)
	disp('Created new table: OperationKeywords');
else
    disp(emsg)
	error('Error creating new table: OperationKeywords')
end

% Cycle through all unique keywords and add them to the OperationKeywords table
K = length(ukws); % the number of unique keywords; the maximum mkw_id index
for k = 1:K
    InsertString = ['INSERT INTO OperationKeywords (Keyword) VALUES (''' ukws{k} ''')'];
    mysql_dbexecute(dbc, InsertString);
end


%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between operation keywords and the operations themselves');

% Create mkwFileRelate Table
CreateString =  SQL_TableCreateString('mkwFileRelate');
[rs,emsg] = mysql_dbexecute(dbc, CreateString);
if ~isempty(rs)
	disp('Created new table: mkwFileRelate');
else
    disp(emsg)
	error('Error creating new table: mkwFileRelate')
end

% Query series table for each keyword
for k = 1:K
    kw = char(ukws{k});
    querystring = ['INSERT INTO mkwFileRelate (mkw_id, m_id) SELECT ' num2str(k) ', m_id FROM Operations ' ...
        	'WHERE (Keywords LIKE ''' kw ',%'' OR Keywords LIKE ''%,' kw ',%'' OR Keywords LIKE ''%,' kw ''' OR' ...
        	' Keywords = ''' kw ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, querystring);
end


%% Go back and write number of occurences to OperationKeywords Table
for k = 1:K
	countstring = ['SELECT COUNT(mkw_id) AS Countme FROM mkwFileRelate WHERE mkw_id = ' num2str(k)];
	[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,countstring);
	if isempty(qrc)
		disp(['Error evaluating COUNT']); disp([emsg]);
	end
	
	% write how many there are back
	updatestring = ['UPDATE OperationKeywords SET NumOccur = ' num2str(qrc{1}) ' WHERE mkw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, updatestring);
end


%% Close database
SQL_closedatabase(dbc)

end