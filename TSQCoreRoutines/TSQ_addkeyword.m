function TSQ_addkeyword(morts,ids,kwadd,dbname,updateall,tostart)
%%% TSQ_addkeyword
% Takes in a keyword and adds them to the ids given
% Ben Fulcher 13/1/2010
% Ben Fulcher 12/5/2010 added updateall option
% Ben Fulcher 7/10/2010 added tostart option

if strcmp(morts,'ts')
	disp('Time Series');
elseif strcmp(morts,'mets')
	disp('Operations');
else
	disp('You must specify either ''ts'' or ''mets''. Exiting.');
	return
end

if nargin < 4
	dbname = '';
end

if nargin < 5 || isempty(updateall)
	updateall = 1; % by default go and update time series keywords tables, etc.
end

if nargin < 6 || isempty(tostart)
    tostart = 0; % add keyword to end
end


%% Open MySQL Database
dbc = SQL_opendatabase(dbname);

id_string = bencat(ids,',');

if strcmp(morts,'ts') % Time Series
	% Get the ts_ids, Keywords from TimeSeries Table
	SelectString = ['SELECT ts_id, FileName, Keywords FROM TimeSeries WHERE ts_id IN (' id_string ')'];
	[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
	if ~isempty(emsg)
		disp('Error selecting from TimeSeries');
		disp(emsg)
		keyboard
	end
	
	tsids = vertcat(qrc{:,1});
	filenames = qrc(:,2);
	keywords = qrc(:,3);
	nget = length(keywords);
	
	% Add kwadd to keywords
	disp(['Adding ' kwadd ' to all entries']);
    if tostart
        for i = 1:nget
            keywords{i} = [kwadd ',' keywords{i}];
        end
    else
        for i = 1:nget
            keywords{i} = [keywords{i} ',' kwadd];
        end
    end
	
	% Display user output
	for i = 1:nget
		fprintf('%s\t%s\n',filenames{i},keywords{i})
	end
	
	% Write back
	reply = input('Sure you want me to write back? Give me an ''n'' if not...','s');
    if strcmp(reply,'n')
        disp('fine -- nothing more will be done on my watch')
        return
    end
	times = zeros(nget,1);
	% do it individually
	for i=1:nget
		tic

		updatestring = ['UPDATE TimeSeries SET ' ...
							'Keywords = ''' keywords{i} ''', ' ...
							'LastModified = NOW() ' ...
							'WHERE ts_id = ' num2str(tsids(i))];
	    [~,emsg] = mysql_dbexecute(dbc, updatestring);
		if ~isempty(emsg)
			disp('Error writing new keywords to TimeSeries');
			disp(emsg)
			keyboard
		end
		
		times(i) = toc;
		if mod(i,floor(nget/5))==0
			disp(['Approximately ' benrighttime(mean(times(1:i))*(nget-i)) ' remaining!']);
		end
	end
	
	% Redo keywords
	if updateall
		disp('Updating TimeSeries Keywords -- create new relation tble tskwRelate');
		SQL_update_tskw(dbname);
	
		% Updating PercentageCalculated ONLY for keywords (should be quick)
		SQL_fillfromresults(tsids,[],[1 1 1],[0 1 0 0],dbname);
	end
	
else % Operations

	% Get the ts_ids, Keywords from TimeSeries Table
	SelectString = ['SELECT m_id, OpName, Keywords FROM Operations WHERE m_id IN (' id_string ')'];
	% makes absolutely certain that ids given match up with keywords obtained
	[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
	if ~isempty(emsg)
		disp('Error selecting from Operations');
		disp(emsg)
		keyboard
	end
	
	mids = vertcat(qrc{:,1});
    opnames = qrc(:,2);
	keywords = qrc(:,3);
	nget = length(keywords);
	
	% Add kwadd to keywords
	disp(['Adding ' kwadd ' to all entries']);
	for i=1:nget
		keywords{i} = [keywords{i} ',' kwadd];
        disp([opnames{i} '   [' keywords{i} ']'])
	end
	
	% Write back
	reply = input('Sure you want me to write back? Give me an ''n'' if not...','s');
    if strcmp(reply,'n')
        disp('fine -- nothing more will be done on my watch')
        return
    end
	times = zeros(nget,1);
	% do it individually
	for i=1:nget
		tic

		updatestring = ['UPDATE Operations SET ' ...
							'Keywords = ''' keywords{i} ''', ' ...
							'LastModified = NOW() ' ...
							'WHERE m_id = ' num2str(mids(i))];
	    [~,emsg] = mysql_dbexecute(dbc, updatestring);
		if ~isempty(emsg)
			disp('Error writing new keywords to TimeSeries');
			disp(emsg)
			keyboard
		end
		
		times(i) = toc;
		if mod(i,floor(nget/5))==0
			disp(['Approximately ' benrighttime(mean(times(1:i))*(nget-i)) ' remaining!']);
		end
	end
	
	% Redo keywords
	if updateall
		disp('Updating Operation Keywords -- create new relation tble mkwRelate');
		SQL_update_mkw(dbname);
	
		% Updating PercentageCalculated, etc. ONLY for keywords (should be quick)
		SQL_fillfromresults([],mids,[1 1 1],[0 0 0 1],dbname);
	end
end

% Close database connection
SQL_closedatabase(dbc)

disp('Done. Not painful at all. I told you to trust me.');

end