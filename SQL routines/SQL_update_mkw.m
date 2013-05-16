%%% SQL_update_mkw
function SQL_update_mkw(dbname)
	
% recreates the keywords tables and their links to metrics
% To be run when operations are either added or removed
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 turned into a function; added dbname input

if nargin<1
	dbname = ''; % opens default specified in SQL_opendatabase
end

%% Open database
dbc = SQL_opendatabase(dbname); % dbc is the database

%% Drop existing tables
% (1) mkwFileRelate
[thetables,qrf,rs,emsg] = mysql_dbquery(dbc,'SHOW TABLES');
if ismember('mkwFileRelate',thetables)
	% alread exists -- drop and recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE mkwFileRelate');
	if ~isempty(rs)
		disp('mkwFileRelate table dropped');
	else
		disp('Error dropping table mkwFileRelate'); disp(emsg);
	end
end

% (2) Operation Keywords
if ismember('operationkeywords',thetables)
	% alread exists -- drop then recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE OperationKeywords');
	if ~isempty(rs)
		disp('Successfully dropped OperationKeywords table');
	else
		disp('Error dropping OperationKeywords table'); disp(emsg);
	end
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
createstring = ['CREATE TABLE OperationKeywords ' ...
				'(mkw_id integer not null auto_increment primary key, Keyword varchar(255), ' ...
					'NumOccur integer, PercentageCalculated float, PercentageGood float, MeanCalcTime float)'];
mysql_dbexecute(dbc, createstring);
disp('Created new table: OperationKeywords');

% Cycle through all unique keywords and add them to the OperationKeywords table
K = length(ukws); % the number of unique keywords; the maximum mkw_id index
for k = 1:K
    insertstring = ['INSERT INTO OperationKeywords (Keyword) VALUES (''' ukws{k} ''')'];
    mysql_dbexecute(dbc, insertstring);
end


%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between operation keywords and the operations themselves');

% Create mkwFileRelate Table
createstring = ['CREATE TABLE mkwFileRelate (mkw_id integer, m_id integer, '  ...
				'FOREIGN KEY (mkw_id) REFERENCES OperationKeywords (mkw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
				'FOREIGN KEY (m_id) REFERENCES Operations (m_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if ~isempty(rs)
	disp('Created new table: mkwFileRelate');
else
	disp('Error creating new table: mkwFileRelate'); disp([emsg]);
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
for k=1:K
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