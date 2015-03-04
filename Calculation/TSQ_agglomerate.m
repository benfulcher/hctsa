% ------------------------------------------------------------------------------
% TSQ_agglomerate
% ------------------------------------------------------------------------------
% 
% Uploads data in the HCTSA_loc file in the current directory back into the
% mySQL database. Should be done to store the result new computations done by
% TSQ_brawn.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TSQ_agglomerate(writeWhat,logToFile,dbname)

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 1
	writeWhat = 'null'; % 'nullerror'
    % find all nulls in the database and write over them if there are values in local files
end
if ~ismember(writeWhat,{'null','error','nullerror'})
    error('Unknown specifier ''%s''',writeWhat)
end
if nargin < 2 || isempty(logToFile)
	logToFile = 0;
end
if nargin < 3
	dbname = '';
end

% ------------------------------------------------------------------------------
%% Open a mySQL database connection
% ------------------------------------------------------------------------------
[dbc, dbname] = SQL_opendatabase(dbname);

% ------------------------------------------------------------------------------
%% Load local files
% ------------------------------------------------------------------------------
%% Read in information from local files
fid = 1; % haha no more logging option...
fprintf(fid,'Loading data from HCTSA_loc.mat...');
load('HCTSA_loc.mat')
fprintf(fid,' Done.\n');

% ------------------------------------------------------------------------------
%% Preliminary definitions
% ------------------------------------------------------------------------------
nts = length(TimeSeries); % Number of time series
nm = length(Operations); % Number of operations
ts_ids_string = BF_cat([TimeSeries.ID],',');
op_ids_string = BF_cat([Operations.ID],',');

% ------------------------------------------------------------------------------
%% Check that nothing has been deleted since the data was loaded...
% ------------------------------------------------------------------------------
% Time series
SelectString = sprintf('SELECT COUNT(ts_id) FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string);
nts_db = mysql_dbquery(dbc,SelectString);
nts_db = nts_db{1};

% Operations
SelectString = sprintf('SELECT COUNT(op_id) FROM Operations WHERE op_id IN (%s)',op_ids_string);
nop_db = mysql_dbquery(dbc,SelectString);
nop_db = nop_db{1};

if (nm == nop_db) && (nts == nts_db)
    fprintf(1,'All local time series and operation ids still exist in %s. This is good.\n',dbname)
else
    if nts_db < nts
        fprintf(1,'There are %u time series that no longer match the database',(nts-nts_db));
    end

    if nop_db < nm
    	fprintf(1,'There are %u operations that no longer match the database',(nm-nop_db));
    end    
    error(['It could be dangerous to write back to a changed database. ' ...
                'You should start a TSQ_prepared from scratch.'])
end

% ------------------------------------------------------------------------------
%% Find the elements that are empty in the database
%       (and hopefully full in the local file)
% ------------------------------------------------------------------------------
% Parts of calculated subsection that are empty in storage
fprintf(1,'Retrieving %s elements from the Results table in %s...',writeWhat,dbname);

switch writeWhat
case 'null'
    % collect nulls in the database
    SelectString = sprintf(['SELECT ts_id, op_id FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND QualityCode IS NULL'],ts_ids_string,op_ids_string);
case 'nullerror'
    % Collect all NULLS and previous errors
    SelectString = sprintf(['SELECT ts_id, op_id, QualityCode FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND (QualityCode IS NULL OR QualityCode = 1)'], ...
        					ts_ids_string,op_ids_string);
case 'error'
    % Collect all previous errors (assume done a TSQ_prepared using 'error' input)
    SelectString = sprintf(['SELECT ts_id, op_id FROM Results WHERE ts_id IN (%s)' ...
    					' AND op_id IN (%s) AND QualityCode = 1'], ...
        					ts_ids_string,op_ids_string);
end

RetrievalTimer = tic; % Time the retrieval (should be fast)
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
if ~isempty(emsg)
    fprintf(1,'\n'); error('Error selecting %s elements from %s',writeWhat,dbname);
elseif isempty(qrc)
    fprintf(1,'\nNo %s elements in this range in the database anymore!\n',writeWhat);
    SQL_closedatabase(dbc); return
else
	fprintf(1,' Retrieved %u entries in %s\n',length(qrc),BF_thetime(toc(RetrievalTimer)));
end
clear RetrievalTimer % Stop timing


ts_id_db = vertcat(qrc{:,1}); % ts_ids (in op_id pairs) of empty database elements in this ts_id/op_id range
op_id_db = vertcat(qrc{:,2}); % op_ids (in ts_id pairs) of empty database elements in this ts_id/op_id range
ndbel = length(ts_id_db);     % Number of database elements to attempt to write back to

% ------------------------------------------------------------------------------
% Give user feedback about what database writing will occur
% ------------------------------------------------------------------------------
switch writeWhat
case 'null'
    fprintf(1,['There are %u NULL entries in Results.\nNow writing calculated ' ...
                    'elements of TS_DataMat to %s...\n'],ndbel,dbname);
case 'error'
    fprintf(1,['There are %u entries in Results (all previous errors) ' ...
                    'that are being written to %s...\n'],ndbel,dbname);
    fprintf(1,['Previous results stored as errors in the database WILL NOT ' ...
                                    'be overwritten with newer errors\n'])
case 'nullerror'
    q_db = qrc(:,3); % empties (NULL) and fatal error (1)
    q_db(cellfun(@isempty,q_db)) = {0}; % turn NULLs to 0s
    q_db = vertcat(q_db{:}); % turn cells to a numeric vector
    % so now nulls in database are labeled 0, and previous errors are labeled 1
    fprintf(1,['There are %u entries in Results (either NULL or previous errors) ' ...
                    'that are being written to %s...\n'],ndbel,dbname);
    fprintf(1,['Note that previous results stored as errors in the database WILL NOT ' ...
                                'be overwritten with newer errors\n'])
    fprintf(1,'However, NULLS will be written over with any result from the local files\n')
end

LocalIndex = zeros(ndbel,2);
LocalIndex(:,1) = arrayfun(@(x)find([TimeSeries.ID] == x,1),ts_id_db); % Indices of rows in local file for each entry in the database
LocalIndex(:,2) = arrayfun(@(x)find([Operations.ID] == x,1),op_id_db); % Indicies of columns in local file for each entry in the database
updateMe = zeros(ndbel,1); % Label iterations that should be written to the database

% ------------------------------------------------------------------------------
% Begin writing each element to the database, one at a time, using UPDATE commands
% ------------------------------------------------------------------------------

writeBackTimer = tic; % Time how long this takes to give user feedback

for i = 1:ndbel	
    
    % Retrieve the elements
    TS_DataMat_ij = TS_DataMat(LocalIndex(i,1),LocalIndex(i,2));
    TS_Quality_ij = TS_Quality(LocalIndex(i,1),LocalIndex(i,2));
    TS_CalcTime_ij = TS_CalcTime(LocalIndex(i,1),LocalIndex(i,2));
    
    switch writeWhat
    case 'null'
        if isfinite(TS_DataMat_ij)
            updateMe(i) = 1; % There is a value in TS_DataMat -- write it back to the NULL entry in the database
        end
    case 'error'
        if isfinite(TS_DataMat_ij) && TS_Quality_ij~=1
            updateMe(i) = 1; % There is a now a non-error value in TS_DataMat previously returned an error (in the database)
        end
    case 'nullerror'
        if isfinite(TS_DataMat_ij) && (q_db(i)==0 || TS_Quality_ij~=1)
            updateMe(i) = 1; % there is a value in TS_DataMat -- write it to the entry in the database
        end
		% (i) Has been calculated and a value stored in TS_DataMat (isfinite()), and 
		% (ii) Either the database entry is NULL or we didn't get an error (prevents writing errors over errors)
    end
	
    if updateMe(i)
        
        if isnan(TS_CalcTime_ij) % happens when there is an error in the code
            TS_CalcTime_string = 'NULL';
        else
            TS_CalcTime_string = sprintf('%f',TS_CalcTime_ij);
        end
            
        % I can't see any way around running lots of single UPDATE commands (for each entry)
    	UpdateString = sprintf(['UPDATE Results SET Output = %19.17g, QualityCode = %u, CalculationTime = %s ' ...
							'WHERE ts_id = %u AND op_id = %u'],TS_DataMat_ij,TS_Quality_ij, ...
							TS_CalcTime_string,ts_id_db(i),op_id_db(i));
        [~,emsg] = mysql_dbexecute(dbc, UpdateString);
        if ~isempty(emsg)
            SQL_closedatabase(dbc) % close the database connection first...
        	error('Error storing (ts_id,op_id) = (%u,%u) to %s??!!\n%s\n', ...
                			[TimeSeries(LocalIndex(i,1)).ID],[Operations(LocalIndex(i,2)).ID],dbname,emsg);
        end
    end

    % Give user feedback on how long is remaining:
    if (i==1)
        fprintf(1,['Based on the first retrieval, this is taking ' ...
                'approximately %s per entry to write to the database.\n'],BF_thetime(toc(writeBackTimer)));
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(toc(writeBackTimer)/i*(ndbel-i)));
    elseif mod(i,floor(ndbel/5))==0 % Give 5 more time updates...
		fprintf(1,['Approximately %s remaining! -- so far %u entries (/ %u possible) have been'  ...
			' written to %s...\n'],BF_thetime(toc(writeBackTimer)/i*(ndbel-i)),sum(updateMe),i,dbname);
	end
end

% ------------------------------------------------------------------------------
% Finished writing back to the database!
% ------------------------------------------------------------------------------

fprintf(1,['So that all seemed to go well -- we wrote %u new calculation results ' ...
                '(/ %u) to the Results table in %s.\n'],sum(updateMe),ndbel,dbname);
fprintf(1,'Writing to the database took at total of %s.\n',BF_thetime(toc(writeBackTimer)));

if any(~updateMe) % Some entries were not written to the database
    fprintf(1,['%u entries were not written (previously-calculated errors) and remain ' ...
                            'awaiting calculation in the database.\n'],sum(~updateMe));
end


clear writeBackTimer % Stop the timer

SQL_closedatabase(dbc) % Close the database connection

% if logToFile
%     fprintf(1,'Logging to file...\n');
%     fn = ['TS_agglomerate_' datestr(now,30) '.log'];
%     fid = fopen(fn,'w','n');
%     disp(['Log file created: ' fn]);
% 
%     fprintf(fid, '%s\n', ['Updated ' num2str(length(tsgoodi)) ' time series: ']);
%     for i=1:length(tsgoodi)
%         fprintf(fid, '%s\n',tsf{tsgoodi(i)});
%     end
% 
%     fprintf(fid, '\n\n\n\n\n%s\n', '******************************');
%     fprintf(fid, '%s\n', ['Updated ' num2str(length(mgoodi)) ' operations: ']);
%     for i=1:length(mgoodi)
%         fprintf(fid, '%s\n',mlab{mgoodi(i)});
%     end
%     fclose(fid);
% end


end