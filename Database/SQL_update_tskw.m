function SQL_update_tskw(dbname)
% Creates keyword associations from comma-delimited columns in the TimeSeries table
% 				   1) splits off a table of unique, seperate keywords
% 				   2) creates association table of keywords to timeseries
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 Added dbname input to specify database

% Note that since table is dropped and recreated, the autonumber columns are guaranteed to be consecutive,
% so we can use the index to know the id (unlike in general, where deletions will change the auto_counter)

if nargin < 1
	dbname = '';
end

%% Open database
dbc = SQL_opendatabase(dbname); % dbc is the database

%% Drop existing Tables
% Will be recreated in what follows
[thetables,qrf,rs,emsg] = mysql_dbquery(dbc,'SHOW TABLES');
% Must be done first because of foreign key constraint
if ismember('tskwFileRelate',thetables)
	% alread exists -- drop and recreate it
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE tskwFileRelate');
	if ~isempty(rs)
		disp('tskwFileRelate dropped');
	else
        disp(emsg)
		error('Error dropping tskwFileRelate')
	end
else
    error('tskwFileRelate doesn''t exist??')
end

if ismember('TimeSeriesKeywords',thetables)
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE TimeSeriesKeywords');
	if ~isempty(rs)
		disp('TimeSeriesKeywords table dropped')
	else
        disp(emsg);
		error('Error dropping TimeSeriesKeywords')
	end
else
    error('TimeSeriesKeywords doesn''t exist??')
end


%% Find all unique keywords strings, split into a table of unique, separate keywords
disp('Looking for unique keywords -- to place in new table TimeSeriesKeywords');
 
SelectString = 'SELECT DISTINCT Keywords FROM TimeSeries';
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

CreateString = SQL_TableCreateString('TimeSeriesKeywords')
[rs,emsg] = mysql_dbexecute(dbc, CreateString);
if ~isempty(rs)
	disp('Created new table: TimeSeriesKeywords')
else
    disp(emsg)
	error('Error creating new table: TimeSeriesKeywords')
end

% Cycle through all unique keywords and add them to the TimeSeriesKeywords table
K = length(ukws); % the number of unique keywords; the maximum tskw_id index
for k = 1:K
    InsertString = ['INSERT INTO TimeSeriesKeywords (Keyword) VALUES (''' ukws{k} ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, InsertString);
	if isempty(rs)
		disp(['Error inserting ' ukws{k}]); keyboard
	end
end

%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between time series keywords and the time series themselves');

CreateString = SQL_TableCreateString('tskwFileRelate');
[rs,emsg] = mysql_dbexecute(dbc, CreateString);
if ~isempty(rs)
	disp('Created new table: tskwFileRelate');
else
	disp('Error creating table: tskwFileRelate'); keyboard
end

% Query series table for each keyword
for k = 1:K
    kw = char(ukws{k});
    InsertString = ['INSERT INTO tskwFileRelate (tskw_id, ts_id) SELECT ' num2str(k) ', ts_id FROM TimeSeries ' ...
        	'WHERE (Keywords LIKE ''' kw ',%'' OR Keywords LIKE ''%,' kw ',%'' OR Keywords LIKE ''%,' kw ''' OR' ...
        	' Keywords = ''' kw ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, InsertString);
	if isempty(rs)
		disp(['Error inserting ' kw]);
	end
end


%% Go back and write number of occurrences to TimeSeriesKeywords Table
for k = 1:K
	UpdateString = ['UPDATE TimeSeriesKeywords SET NumOccur = ' ...
							'(SELECT COUNT(tskw_id) AS Countme FROM tskwFileRelate WHERE tskw_id = ' num2str(k) ') ' ...
					'WHERE tskw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);
end

%% Populate the MeanLength column in TimeSeriesKeywords
% Mean length of time series with this keyword
% Requires that lengths are already included as metadata in the TimeSeries table
for k = 1:K
	UpdateString = ['UPDATE TimeSeriesKeywords SET MeanLength = ' ...
						'(SELECT AVG(Length) FROM TimeSeries WHERE ts_id IN (' ...
											'SELECT ts_id FROM tskwFileRelate WHERE tskw_id = ' num2str(k) ')) ' ...
					'WHERE tskw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);						
	if ~isempty(emsg)
		disp(['Error updating TimeSeriesKeywords with MeanLength']); keyboard
	end
end


%% Close database
SQL_closedatabase(dbc)

end