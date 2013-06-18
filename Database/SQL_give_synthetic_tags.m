% Quick script to populate IsSynthetic column from keywords
% Ben Fulcher, May 2013

dbc = SQL_opendatabase;
SelectString = 'SELECT ts_id FROM TimeSeries WHERE Keywords LIKE ''%synthetic%''';
[tsids,~,~,emsg] = mysql_dbquery(dbc,SelectString);
alltsids = mysql_dbquery(dbc,'SELECT ts_id FROM TimeSeries')
synthetictsids = vertcat(tsids{:});
alltsids = vertcat(alltsids{:});
nts = length(alltsids);

for i = 1:nts
    if ismember(alltsids(i),synthetictsids)
        UpdateString = ['UPDATE TimeSeries SET IsSynthetic = 1 WHERE ts_id = ' num2str(alltsids(i))];
    else
        UpdateString = ['UPDATE TimeSeries SET IsSynthetic = 0 WHERE ts_id = ' num2str(alltsids(i))];
    end
	[rs, emsg] = mysql_dbexecute(dbc, UpdateString);
    if ~isempty(emsg)
        disp(['Error updating ' num2str(alltsids(i))])
    end
end
SQL_closedatabase(dbc)
disp('Database closed, IsSynthetic column in the TimeSeries table updated ///');