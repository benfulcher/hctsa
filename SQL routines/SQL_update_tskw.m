%%% SQL_update_tskw
function SQL_update_tskw(dbname)
	
% creates keyword associations from comma-delimited columns in the TimeSeries table
% 				   1) splits off a table of unique, seperate keywords
% 				   2) creates association table of keywords to timeseries
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 Added dbname input to specify database

% Note that since table is dropped and recreated, the autonumber columns are guaranteed to be consecutive,
% so we can use the index to know the id (unlike in general, where deletions will change the auto_counter)

if nargin<1
	dbname = '';
end

%% Open database
dbc = SQL_opendatabase(dbname); % dbc is the database


%% Drop existing Tables
% Will be recreated in what follows
[thetables,qrf,rs,emsg] = mysql_dbquery(dbc,'show tables');
% Must be done first because of foreign key constraint
if ismember('tskwFileRelate',thetables)
	% alread exists -- drop and recreate it
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE tskwFileRelate');
	if ~isempty(rs)
		disp('tskwFileRelate dropped');
	else
		disp(['Error dropping tskwFileRelate']); keyboard
	end
end
if ismember('timeserieskeywords',thetables)
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE TimeSeriesKeywords');
	if ~isempty(rs)
		disp('TimeSeriesKeywords table dropped');
	else
		disp(['Error dropping TimeSeriesKeywords']); keyboard
	end
end


%% Find all unique keywords strings, split into a table of unique, separate keywords
disp('Looking for unique keywords -- to place in new table TimeSeriesKeywords');
 
selectstring = 'SELECT DISTINCT Keywords FROM TimeSeries';
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


createstring = ['CREATE TABLE TimeSeriesKeywords (tskw_id integer auto_increment primary key, Keyword varchar(50), ' ...
						'NumOccur integer, PercentageCalculated float, PercentageGood float, MeanCalcTime float, MeanLength integer)'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if ~isempty(rs)
	disp('Created new table: TimeSeriesKeywords');
else
	disp('Error creating new table: TimeSeriesKeywords'); keyboard
end

% Cycle through all unique keywords and add them to the TimeSeriesKeywords table
K = length(ukws); % the number of unique keywords; the maximum tskw_id index
for k = 1:K
    insertstring = ['INSERT INTO TimeSeriesKeywords (Keyword) VALUES (''' ukws{k} ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
	if isempty(rs)
		disp(['Error inserting ' ukws{k}]); keyboard
	end
end

%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between time series keywords and the time series themselves');

createstring = ['CREATE TABLE tskwFileRelate (tskw_id integer, ts_id integer, ' ...
				'FOREIGN KEY (tskw_id) REFERENCES TimeSeriesKeywords(tskw_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
				'FOREIGN KEY (ts_id) REFERENCES TimeSeries(ts_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if ~isempty(rs)
	disp('Created new table: tskwFileRelate');
else
	disp('Error creating table: tskwFileRelate'); keyboard
end

% Query series table for each keyword
for k = 1:K
    kw = char(ukws{k});
    insertstring = ['INSERT INTO tskwFileRelate (tskw_id, ts_id) SELECT ' num2str(k) ', ts_id FROM TimeSeries ' ...
        	'WHERE (Keywords LIKE ''' kw ',%'' OR Keywords LIKE ''%,' kw ',%'' OR Keywords LIKE ''%,' kw ''' OR' ...
        	' Keywords = ''' kw ''')'];
    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
	if isempty(rs)
		disp(['Error inserting ' kw]);
	end
end


%% Go back and write number of occurrences to TimeSeriesKeywords Table
for k=1:K
	% countstring = ['SELECT COUNT(tskw_id) AS Countme FROM tskwFileRelate WHERE tskw_id = ' num2str(k)];
	% [qrc,qrf,rs,emsg] = mysql_dbquery(dbc,countstring);	
	% % write how many there are back
	% updatestring = ['UPDATE TimeSeriesKeywords SET NumOccur = ' num2str(qrc{1}) ' WHERE tskw_id = ' num2str(k)];
	% [rs,emsg] = mysql_dbexecute(dbc, updatestring);

	updatestring = ['UPDATE TimeSeriesKeywords SET NumOccur = ' ...
							'(SELECT COUNT(tskw_id) AS Countme FROM tskwFileRelate WHERE tskw_id = ' num2str(k) ') ' ...
					'WHERE tskw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, updatestring);
end

%% Fill Mean Length column
% Mean length of time series with this keyword
% Requires that lengths are already included as metadata in the TimeSeries table
for k=1:K
	% selectstring = ['SELECT AVG(Length) FROM TimeSeries WHERE ts_id IN (' ...
						% 'SELECT ts_id FROM tskwFileRelate WHERE tskw_id = ' num2str(k) ')'];
	updatestring = ['UPDATE TimeSeriesKeywords SET MeanLength = ' ...
						'(SELECT AVG(Length) FROM TimeSeries WHERE ts_id IN (' ...
											'SELECT ts_id FROM tskwFileRelate WHERE tskw_id = ' num2str(k) ')) ' ...
					'WHERE tskw_id = ' num2str(k)];
	[rs,emsg] = mysql_dbexecute(dbc, updatestring);						
	if ~isempty(emsg)
		disp(['Error updating TimeSeriesKeywords with MeanLength']); keyboard
	end
end



%% Close database
SQL_closedatabase(dbc)

end