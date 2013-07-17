%%% SQL_initiate_results
% Initiates the results matrix, given that the Operations and TimeSeries databases are already up and running
% Ben Fulcher 29/11/09

%% Open database
[dbc dbname] = SQL_opendatabase; % the default database, called dbname, is opened as dbc

%% Check Results table doesn't already exist -- if it does ask to drop it
thetables = mysql_dbquery(dbc,'Show Tables');

% Check TimeSeries and Operations tables exist ok (Results table uses foreign keys from TimeSeries and Operations tables)
if ~ismember('Operations',thetables) || ~ismember('TimeSeries',thetables)
	disp(['You''re going to need the TimeSeries and Operations tables if you want me to populate the ' ...
				'Results table in ' dbname])
	SQL_closedatabase(dbc)
	beep; return
end
% If Results Table exists, drop it
if ismember('Results',thetables)
	input(['Do you really want me to DROP the ''Results'' table from ' dbname '??! ' ...
				'This will delete everything currently in the Results table. Reply with ''y'' to do this...?'],'s');
	if ~strcmp(reply,'y')
		disp('Ok, better to be safe than sorry.')
		SQL_closedatabase(dbc)
		% get out now if you don't want want to trash your existing Results table
		beep; return
	end
	mysql_dbexecute(dbc, 'DROP TABLE Results')
end
%% Create Results table
% Columns:
% 1) ts_id (unique integer identified from TimeSeries)
% 2) m_id (unique integer identified from Operations)
% 3) Output (double)
% 4) QualityCode (positive integer)
% 5) CalculationTime (float)
% 6) LastModified (datetime)
% [[could add a Machine label for machine calculated on]]
createstring = ['CREATE TABLE Results (ts_id integer, m_id integer, Output double, QualityCode integer unsigned, CalculationTime float, LastModified datetime, ' ...
				'FOREIGN KEY (ts_id) REFERENCES TimeSeries (ts_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
				'FOREIGN KEY (m_id) REFERENCES Operations (m_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if isempty(emsg)
	disp(['Results Table successfully initiated into ' dbname]);
else
	disp(['Error setting up Results Table in ' dbname]); keyboard
end

%% Populate Results table with all combinations of time series and operations to Results table
% Relevant user information/warnings:
disp(['Currently adding all combinations of time series and operations as blank (NULL) entries' ...
			'\n in the Results Table of ' dbname]);
disp(['Note that if the TimeSeries and Operations tables are already populated, this could take a long time.\n Do not interrupt.']);
% Count number of Time series and Operations in the database:
[nm,qrf,rs,emsg] = mysql_dbquery(dbc,'SELECT COUNT(m_id) as nm FROM Operations'); nm = nm{1};
[nts,qrf,rs,emsg] = mysql_dbquery(dbc,'SELECT COUNT(ts_id) as nts FROM TimeSeries'); nts = nts{1};
disp(['Ok, it looks like we have around ' num2str(nm) ' operations and ' num2str(nts) ' time series\n on our hands... Let''s go :)']);

%% -- Do it using the mySQL JOIN command
tic
insertstring = 'INSERT INTO Results (ts_id, m_id) SELECT ts_id, m_id FROM TimeSeries CROSS JOIN Operations';
[rs,emsg] = mysql_dbexecute(dbc, insertstring);
if isempty(rs)
	disp(['Error initializing Results table into ' dbname]);
else
	disp(['Results table successfully initialized into ' dbname]);
end

disp(['DONE!! Successfully!! Really?! Wow!! I bet that took ages. In fact it only took ' BF_thetime(toc)]);

%% -- MATLAB WAY (one-at-a-time)
% SelectString = ['SELECT ts_id from TimeSeries'];
% [all_ts_ids,qrf,rs,emsg] = mysql_dbquery(dbc, SelectString);
% all_ts_ids = vertcat(all_ts_ids{:});
% nts = length(all_ts_ids);
% % all ts_ids in TimeSeries are in the column vector all_ts_ids (ascending integers)
% 
% SelectString = ['SELECT m_id from Operations'];
% [all_m_ids,qrf,rs,emsg] = mysql_dbquery(dbc, SelectString);
% all_m_ids = vertcat(all_m_ids{:});
% nm = length(all_m_ids);
% % all m_ids in Operations are in all_m_ids (ascending integers)
% 
% for i=1:nts
% 	tic
% 	for j=1:nm
% 		insertstring = ['INSERT INTO Results (ts_id, m_id) VALUES (' num2str(i) ',' num2str(j) ')'];
% 	    mysql_dbexecute(dbc, insertstring);
% 	end
% 	disp([num2str(i) ' || ' BF_thetime(toc*(nts-i))]);
% end

%% Close database
SQL_closedatabase(dbc)