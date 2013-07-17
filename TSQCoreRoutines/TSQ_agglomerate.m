function TSQ_agglomerate(writewhat,tolog,dbname)
% agglomerates the TS_loc file in the current directory
% into the database
% Needs files TS_loc, TS_guide_ts, TS_guide_mets in current directory
% 12/1/2010 ~Ben Fulcher~ Complete rewrite -- checks what parts of the storage are unwritten, and writes
% 			only to that section. So it doesn't matter if calculated bits are inconsistent, or if some 
% 			bits have subsequently been calculated by other machines in parallel -- it will only write to
% 			empty parts of the storage which have been calculated here. This makes it quicker (only retrieves part
% 			of storage that is empty)
% 			Also added default logging to file


%% Check Inputs
if nargin < 1
	writewhat = 'null'; % 'nullerror'
    % find all nulls in the database and write over them if there are values in local files
end
if ~ismember(writewhat,{'null','nullerror'})
    error(['Unknown specifier ''' writewhat ''''])
end
if nargin < 2 || isempty(tolog)
	tolog = 0;
end
if nargin < 3
	dbname = '';
end

%% Open a mySQL database connection
[dbc,dbname] = SQL_opendatabase(dbname);

%% Load local files
fprintf(1,'Loading local files...');
load TS_loc.mat TS_loc; fprintf(1,' TS_loc');
load TS_loc_ct.mat TS_loc_ct; fprintf(1,', TS_loc_ct');
load TS_loc_q.mat TS_loc_q; fprintf(1,', TS_loc_q');
load TS_loc_guides.mat ts_ids_keep tsf tskw tsl m_ids_keep mcode ...
            mlab mkw mlink Mmid Mmlab Mmcode % mpoint % maybe don't need to load all of these...?
fprintf(1,', TS_loc_guides. All loaded.\n');

%% Preliminary definitions
nts = length(ts_ids_keep); % number of time series
nm = length(m_ids_keep); % number of operations
ts_ids_keep_string = BF_cat(ts_ids_keep,',');
m_ids_keep_string = BF_cat(m_ids_keep,',');

%% Check that nothing has been deleted in the meantime...
% time series
SelectString = sprintf('SELECT COUNT(ts_id) FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_keep_string);
nts_db = mysql_dbquery(dbc,SelectString);
nts_db = nts_db{1};

% Operations
SelectString = sprintf('SELECT COUNT(m_id) FROM Operations WHERE m_id IN (%s)',m_ids_keep_string);
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


%% Find the elements that are empty in storage (and hopefully full in the local file)
% Parts of calculated subsection that are empty in storage
fprintf(1,'Retrieving %s elements from the Results table in %s...',writewhat,dbname);

switch writewhat
case 'null'
    % collect nulls in the database
    SelectString = sprintf(['SELECT ts_id, m_id FROM Results WHERE ts_id IN (%s)' ...
    					' AND m_id IN (%s) AND QualityCode IS NULL'],ts_ids_keep_string,m_ids_keep_string);
case 'nullerror'
    % collect all NULLS and previous errors
    SelectString = sprintf(['SELECT ts_id, m_id, QualityCode FROM Results WHERE ts_id IN (%s)' ...
    					' AND m_id IN (%s) AND (QualityCode IS NULL OR QualityCode = 1)'], ...
        					ts_ids_keep_string,m_ids_keep_string);
end
tic
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
if ~isempty(emsg)
    fprintf(1,'\n'); error('Error selecting %s elements from %s',writewhat,dbname);
elseif isempty(qrc)
    fprintf(1,'\nNo %s elements in this range in the database anymore!\n',writewhat);
    return
else
	fprintf(1,' Retrieved %u entries in %s\n',length(qrc),BF_thetime(toc));
end

ts_id_db = vertcat(qrc{:,1}); % ts_ids (in m_id pairs) of empty database elements in this ts_id/m_id range
m_id_db = vertcat(qrc{:,2}); % m_ids (in ts_id pairs) of empty database elements in this ts_id/m_id range
ndbel = length(ts_id_db); % number of database elements to (maybe) write back to

switch writewhat
case 'null'
    fprintf(1,['There are %u NULL entries in Results.\nWill now attempt to write calculated elements of TS_loc ' ...
                    'into these elements of %s...\n'],ndbel,dbname);
case 'nullerror'
    q_db = qrc(:,3); % empties (NULL) and fatal error (1)
    q_db(cellfun(@isempty,q_db)) = {0}; % turn NULLs to 0s
    q_db = vertcat(q_db{:}); % turn cells to a numeric vector
    % so now nulls in database are labeled 0, and previous errors are labeled 1
    fprintf(1,['There are %u entries in Results (either NULL or previous errors) ' ...
                    'that are being written to %s...\n'],ndbel,dbname);
    fprintf(1,'Note that previous results stored as errors in the database will not be overwritten with newer errors\n')
    fprintf(1,'However, NULLS will be written over with any result from the local files\n')
end

times = zeros(ndbel,1); % time each iteration
loci = zeros(ndbel,2);
loci(:,1) = arrayfun(@(x)find(ts_ids_keep == x,1),ts_id_db); % indices of rows in local file for each entry in the database
loci(:,2) = arrayfun(@(x)find(m_ids_keep == x,1),m_id_db); % indicies of columns in local file for each entry in the database
updated = zeros(ndbel,1); % label when an iteration writes successfully to the database
for i = 1:ndbel
	tic
    
    % retrieve the elements
    TS_loc_ij = TS_loc(loci(i,1),loci(i,2));
    TS_loc_q_ij = TS_loc_q(loci(i,1),loci(i,2));
    TS_loc_ct_ij = TS_loc_ct(loci(i,1),loci(i,2));
    
    switch writewhat
    case 'null'
        if isfinite(TS_loc_ij)
            updated(i) = 1; % there is a value in TS_loc -- write it back to the NULL entry in the database
        end
    case 'nullerror'
        if isfinite(TS_loc_ij) && (q_db(i)==0 || TS_loc_q_ij~=1)
            updated(i) = 1;
        end
		% (i) Has been calculated and a value stored in TS_loc (isfinite()), and 
		% (ii) either the database entry is NULL or we didn't get an error (prevents writing errors over errors)
    end
	
    if updated(i)
        
        if isnan(TS_loc_ct_ij) % happens when there is an error in the code
            TS_loc_ct_string = 'NULL';
        else
            TS_loc_ct_string = sprintf('%f',TS_loc_ct_ij);
        end
            
        % I can't see any way around running lots of single UPDATE commands (for each entry)
    	UpdateString = sprintf(['UPDATE Results SET Output = %19.17g, QualityCode = %u, CalculationTime = %s ' ...
        							'WHERE ts_id = %u AND m_id = %u'],TS_loc_ij,TS_loc_q_ij, ...
            							TS_loc_ct_string,ts_id_db(i),m_id_db(i));
        [~,emsg] = mysql_dbexecute(dbc, UpdateString);
        if ~isempty(emsg)
        	fprintf(1,'\nError storing (ts_id,m_id) = (%u,%u) to %s??!!', ...
                			ts_ids_keep(loci(i,1)),m_ids_keep(loci(i,2)),dbname);
            fprintf(1,'%s\n',emsg);
            keyboard
        end
    end

	times(i) = toc;
	if mod(i,floor(ndbel/5))==0
		fprintf(1,['Approximately %s remaining! -- so far %u entries'  ...
			' (/ %u possible) have been written to %s...\n'],BF_thetime(mean(times(1:i))*(ndbel-i)),sum(updated),i,dbname);
	end
end

fprintf(1,'Well that seemed to go ok -- we wrote %u new calculation results (/ %u) to the Results table in %s\n',sum(updated),nts*nm,dbname);
if any(~updated) % some were not written
	fprintf(1,'%u entries were not written (recurring fatal errors) and remain awaiting calculation in the database\n',sum(~updated));
end
SQL_closedatabase(dbc) % close database connection

% if tolog
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