% SQL_update_opkw
% 
% Recreates the keywords tables and their links to operations
% To be run when operations are either added or removed
% 
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 turned into a function; added dbname input
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Ben D. Fulcher <ben.d.fulcher@gmail.com>
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function SQL_update_opkw(dbname)

if nargin < 1
	dbname = ''; % opens default specified in SQL_opendatabase
end

%% Open database
[dbc,dbname] = SQL_opendatabase(dbname); % dbc is the database

%% Drop existing tables
% (1) OpKeywordsRelate
[thetables,~,~,emsg] = mysql_dbquery(dbc,'SHOW TABLES');
if any(~isempty(regexp('opkeywordsrelate',thetables,'ignorecase')))
	% alread exists -- drop and recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE OpKeywordsRelate');
	if ~isempty(rs)
		disp('OpKeywordsRelate table dropped');
	else
        disp(emsg);
		error('Error dropping table OpKeywordsRelate');
	end
else
    fprintf(1,'No OpKeywordsRelate table in %s\n',dbname);
end

% (2) Operation Keywords
if any(~isempty(regexp('OperationKeywords',thetables,'ignorecase')))
	% alread exists -- drop then recreate
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE OperationKeywords');
	if ~isempty(rs)
		disp('Successfully dropped OperationKeywords table');
	else
        disp(emsg);
		error('Error dropping OperationKeywords table');
	end
else
    fprintf(1,'No OperationKeywords table in %s\n',dbname);
end

%% Find all unique keywords strings, split into a table of unique, separate keywords
disp('Looking for unique keywords -- to place in new table OperationKeywords');
 
SelectString = 'SELECT DISTINCT Keywords FROM Operations';
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);

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
K = length(ukws); % the number of unique keywords; the maximum opkw_id index
for k = 1:K
    InsertString = ['INSERT INTO OperationKeywords (Keyword) VALUES (''' ukws{k} ''')'];
    mysql_dbexecute(dbc, InsertString);
end


%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between operation keywords and the operations themselves');

% Create OpKeywordsRelate Table
CreateString =  SQL_TableCreateString('OpKeywordsRelate');
[rs,emsg] = mysql_dbexecute(dbc, CreateString);
if ~isempty(rs)
	disp('Created new table: OpKeywordsRelate');
else
    disp(emsg)
	error('Error creating new table: OpKeywordsRelate')
end

% Query series table for each keyword
for k = 1:K
    kw = char(ukws{k});
    querystring = ['INSERT INTO OpKeywordsRelate (opkw_id, op_id) SELECT ' num2str(k) ', op_id FROM Operations ' ...
        	'WHERE (Keywords LIKE ''' kw ',%'' OR Keywords LIKE ''%,' kw ',%'' OR Keywords LIKE ''%,' kw ''' OR' ...
        	' Keywords = ''' kw ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, querystring);
end


%% Go back and write number of occurences to OperationKeywords Table
for k = 1:K
	countstring = ['SELECT COUNT(opkw_id) AS Countme FROM OpKeywordsRelate WHERE opkw_id = ' num2str(k)];
	[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,countstring);
	if isempty(qrc)
		disp(['Error evaluating COUNT']); disp([emsg]);
	end
	
	% write how many there are back
	updatestring = ['UPDATE OperationKeywords SET NumOccur = ' num2str(qrc{1}) ' WHERE opkw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, updatestring);
end


%% Close database
SQL_closedatabase(dbc)

end