% function SQL_update_tskw(dbname)
% Creates keyword associations from comma-delimited columns in the TimeSeries table
% Ben Fulcher 24/11/09 (based on code inherited from Max Little)
% Ben Fulcher 11/1/10 Added dbname input to specify database
% Ben Fulcher June 2013. Rewritten to assume just adding incrementally to database -- just adds the new entries

% Note that since table is empty, the autonumber columns are guaranteed to be consecutive so we can use the index to know the id (unlike in general, where deletions will change the auto_counter)


%% Open default database
dbc = SQL_opendatabase; % dbc is the database


%% Find all unique keywords strings, split into a table of unique, separate keywords
disp('Looking for unique keywords in the TimeSeries Table to place in table TimeSeriesKeywords');

SelectString = 'SELECT DISTINCT Keywords FROM TimeSeries';
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);

ukws = {};
for i = 1:length(qrc)
    kws = regexp(qrc{i},',','split','ignorecase');
    for j = 1:length(kws)
        if ~ismember(kws{j},ukws) % add it to ukws
            ukws{end+1} = kws{j};
        end
    end
end

nkw = length(ukws); % the number of unique keywords; the maximum tskw_id index
fprintf(1,'I just found %g unique keywords\n',nkw)

% ukws = unique(splitkws); % cell of unique keyword strings
%% Drop and recreate time series keyword tables
% [thetables,qrf,rs,emsg] = mysql_dbquery(dbc,'SHOW TABLES');
% % Must be done first because of foreign key constraint
% if ismember('TsKeywordsRelate',thetables)
%     % already exists -- truncate it
%     [rs,emsg] = mysql_dbexecute(dbc, 'TRUNCATE TABLE TsKeywordsRelate');
%     if ~isempty(rs)
%         disp('TsKeywordsRelate Table truncated');
%     else
%         disp(emsg)
%         error('Error truncating TsKeywordsRelate')
%     end
%     % % Recreate it
%     % CreateString = SQL_TableCreateString('TsKeywordsRelate');
%     % [rs,emsg] = mysql_dbexecute(dbc, CreateString);
%     % if ~isempty(rs)
%     %     disp('Created new table: TsKeywordsRelate');
%     % else
%     %     disp('Error creating table: TsKeywordsRelate'); keyboard
%     % end
% else
%     error('The TsKeywordsRelate table doesn''t exist??')
% end
% 
% if ismember('TimeSeriesKeywords',thetables)
%     [rs,emsg] = mysql_dbexecute(dbc, 'TRUNCATE TABLE TimeSeriesKeywords');
%     if ~isempty(rs)
%         disp('TimeSeriesKeywords table truncated')
%     else
%         disp(emsg);
%         error('Error truncating TimeSeriesKeywords')
%     end
%     % % Recreate it
%     % CreateString = SQL_TableCreateString('TimeSeriesKeywords');
%     % [rs,emsg] = mysql_dbexecute(dbc, CreateString);
%     % if ~isempty(rs)
%     %     disp('Created new table: TimeSeriesKeywords')
%     % else
%     %     disp(emsg)
%     %     error('Error creating new table: TimeSeriesKeywords')
%     % end
% else
%     error('The TimeSeriesKeywords table doesn''t exist??')
% end


% Cycle through all unique keywords and add them to the TimeSeriesKeywords table
for k = 1:nkw
    InsertString = sprintf('INSERT INTO TimeSeriesKeywords (Keyword) VALUES (''%s'')',ukws{k});
    [rs,emsg] = mysql_dbexecute(dbc, InsertString);
	if ~isempty(emsg)
		fprintf(1,'Error inserting %s\n',ukws{k});
        fprintf(1,'%s\n',emsg)
        keyboard
	end
end
fprintf(1,'Just filled TimeSeriesKeywords Table with unique keywords\n')

%% Associate primary keys of keywords and series
disp('Now creating and filling the association table between time series keywords and the time series themselves');

% Get tskw_ids
SelectString = 'SELECT tskw_id FROM TimeSeriesKeywords';
tskw_ids = mysql_dbquery(dbc,SelectString);
tskw_ids = vertcat(tskw_ids{:});

% Query series table for each keyword
for k = 1:nkw
    InsertString = sprintf(['INSERT INTO TsKeywordsRelate (tskw_id, ts_id) SELECT %u, ts_id FROM TimeSeries ' ...
        	'WHERE (Keywords LIKE ''%s,%%'' OR Keywords LIKE ''%%,%s,%%'' OR Keywords LIKE ''%%,%s'' OR' ...
        	' Keywords = ''%s'')'],tskw_ids(k),ukws{k},ukws{k},ukws{k},ukws{k});
    [rs,emsg] = mysql_dbexecute(dbc, InsertString);
	if ~isempty(emsg)
		fprintf(1,'Error inserting %s\n',ukws{k});
        fprintf(1,'%s\n',emsg)
        keyboard
	end
end
SelectString = 'SELECT COUNT(*) FROM TsKeywordsRelate';
nentries = mysql_dbquery(dbc,SelectString);
fprintf(1,'Just inserted %u keyword relations for %u Keywords into the TsKeywordsRelate table\n',nentries{1},nkw)

%% Go back and write number of occurrences to TimeSeriesKeywords Table
for k = 1:nkw
	UpdateString = sprintf(['UPDATE TimeSeriesKeywords SET NumOccur = ' ...
						'(SELECT COUNT(*) FROM TsKeywordsRelate WHERE tskw_id = %u) ' ...
    					'WHERE tskw_id = %u'],tskw_ids(k),tskw_ids(k));
	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);
end
fprintf(1,'Just updated the number of occurances of the %u keywords into the TsKeywordsRelate table\n',nkw)


%% Populate the MeanLength column in TimeSeriesKeywords
% Mean length of time series with this keyword
% Requires that lengths are already included as metadata in the TimeSeries table
for k = 1:nkw
	UpdateString = sprintf(['UPDATE TimeSeriesKeywords SET MeanLength = ' ...
						'(SELECT AVG(Length) FROM TimeSeries WHERE ts_id IN (' ...
								'SELECT ts_id FROM TsKeywordsRelate WHERE tskw_id = %u)) ' ...
    					'WHERE tskw_id = %u'],tskw_ids(k),tskw_ids(k));
	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);						
	if ~isempty(emsg)
		fprintf(1,'Error updating TimeSeriesKeywords with MeanLength\n%s\n',emsg);
        keyboard
	end
end
fprintf(1,'Just updated the mean length of time series for each of the %u keywords into the TsKeywordsRelate table\n',nkw)


%% Close database
SQL_closedatabase(dbc)
