function TSQ_agglomerate(tolog,dbname)
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
if nargin<1 || isempty(tolog)
	tolog = 0;
end
if nargin<2
	dbname = '';
end

%% Open MySQL Database
dbc = SQL_opendatabase(dbname);

%% Load files
% Load in necessary files for both reading and writing
% For reading: LOCAL FILES

disp('Loading TS_loc');
load TS_loc.mat TS_loc

disp('Loading TS_loc_q');
load TS_loc_q.mat TS_loc_q

disp('Loading TS_loc_ct');
load TS_loc_ct.mat TS_loc_ct

% All of this probably unnecessary:
disp('Loading TS_loc_guides');
load TS_loc_guides.mat ts_ids_keep tsf tskw tsl m_ids_keep mcode mlab mkw mpoint mlink Mmid Mmlab Mmcode

%% Preliminary definitions
nts = length(tsf);
nm = length(mlab);

%% Check that labels are still good, i.e., index structure is maintained...
% First check that tsf obtained from the database is the same as the tsf in the guide
% This should always be the case, given the autonumbering system used in the mySQL database
% Actually, some may have been deleted in the meantime
disp(['Checking time series still match up']);
ts_ids_keep_string = bencat(ts_ids_keep,',');
selectstring = ['SELECT FileName FROM TimeSeries WHERE ts_id IN (' ts_ids_keep_string ')'];
[tsf_check,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
matches = strcmp(tsf,tsf_check);
tsgoodi = find(matches);
ntsgood = length(tsgoodi);
if ntsgood<nts
	disp(['There are ' num2str(nts-ntsgood) ' time series that no longer match the database']);
end

% Check that mlab still matches up with database
disp(['Checking that metrics still match up']);
m_ids_keep_string = bencat(m_ids_keep,',');
selectstring = ['SELECT OpName FROM Operations WHERE m_id IN (' m_ids_keep_string ')'];
[mlab_check,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
matches = strcmp(mlab,mlab_check);
mgoodi = find(matches);
nmgood = length(mgoodi);
if nmgood<nm
	disp(['There are ' num2str(nm-nmgood) ' metrics that no longer match the database']);
end

if nm==nmgood && nts==ntsgood
    disp('All local time series and metrics match up with the storage file indicies. This is a good thing')
else
    if nm<nmgood
        disp('SAVE ONLY THE GOOD OPERATIONS?')
		% i.e. use TS_loc(:,mgoodi) 
    end
    if nts<ntsgood
        disp('SAVE ONLY THE GOOD TIME SERIES?')
		% i.e. use TS_loc(tsgoodi,:)
    end
	disp('Nope, I can''t take this right now. I''m leaving!');
	return
    input('Control-C now or forever hold your peace...')
end


% %% Check no specials in any matricies
% We check for this now
% if any(~isfinite(TS_loc))
% 	disp(['WHAT ARE NANS/INFS DOING IN TS_LOC??!!?']); return
% end
% if any(~isfinite(TS_loc_q)) % this is ok -- we check for this
% 	disp(['WHAT ARE NANS/INFS DOING IN TS_LOC_Q??!!?']); return
% end
% if any(~isfinite(TS_loc_ct))
% 	disp(['WHAT ARE NANS/INFS DOING IN TS_LOC_CT??!!?']); return
% end


%% Find the elements that are empty in storage but full in the local file
% Parts of calculated subsection that are empty in storage
disp(['Retrieving empty bits of storage in the given range... Be patient... We''re timing it, if that helps...']);
tic
selectstring = ['SELECT ts_id, m_id, Quality FROM Results WHERE ts_id IN (' ts_ids_keep_string ')' ...
					' AND m_id IN (' m_ids_keep_string ' ) AND (Quality IS NULL OR Quality = 1)'];
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
if ~isempty(emsg)
	disp(['Error selecting empty elements from the database']);
	disp(emsg); keyboard
elseif isempty(qrc)
	disp(['No empty elements in this range in the database anymore! There''s a surprise!! Exiting.']); return
else
	disp(['Retrieved!! Took ' benrighttime(toc)]);
end

ts_id_empties = vertcat(qrc{:,1}); % ts_ids (in m_id pairs) of empty database elements in this ts_id/m_id range
m_id_empties = vertcat(qrc{:,2}); % m_ids (in ts_id pairs) of empty database elements in this ts_id/m_id range
nempties = length(ts_id_empties);

qualities = qrc(:,3); % empties (null) and 1s (fatal errors)
qualities(cellfun(@isempty,qualities)) = {NaN}; % turn NULLs in NaNs
qualities = vertcat(qualities{:});

% % Convert empty entries to NaNs
% nothings = find(cellfun(@isempty,qualities)); % nothings: where empty entries
% if ~isempty(nothings)
% 	qualitites(nothings) = {NaN};
% end

disp(['There are ' num2str(nempties) ' entries in the database in this range that await filling... Let''s do it...!']);

times = zeros(nempties,1);
updatecounter = 0;
for i = 1:nempties
	tic
	ind_i = find(ts_ids_keep == ts_id_empties(i)); % the row index in TS_loc_ / ts_ids_keep
	ind_j = find(m_ids_keep == m_id_empties(i)); % the column index in TS_loc_ / m_ids_keep
	% could maybe do this all at once with an arrayfun before the loop?
	
	if isfinite(TS_loc_q(ind_i,ind_j)) && isfinite(TS_loc(ind_i,ind_j)) && (isnan(qualities(i)) || TS_loc_q(ind_i,ind_j)~=1)
		% Has to have been attempted to be calculated (isfinite()); and TS_loc hasn't mistakenly been stored with a NaN/Inf
		% (it never should be) and will then write back if either
		% the store has nothing (isnan()) or if we have a non-fatal error to write back. i.e., will NOT write
		% back over a fatal error that is again fatal.
		updatecounter = updatecounter + 1;
		
		% Calcultion time -- store as NULL if NaN
		if isnan(TS_loc_ct(ind_i,ind_j))
			string_ct = 'NULL';
		else
			string_ct = num2str(TS_loc_ct(ind_i,ind_j));
		end

		updatestring = ['UPDATE Results SET ' ...
							'Output = ' num2str(TS_loc(ind_i,ind_j),'%19.17g') ', ' ...
							'Quality = ' num2str(TS_loc_q(ind_i,ind_j)) ', ' ...
							'CalculationTime = ' string_ct ', ' ...
							'LastModified = NOW() ' ...
							'WHERE ts_id = ' num2str(ts_ids_keep(ind_i)) ' AND m_id = ' num2str(m_ids_keep(ind_j))];

	    [rs,emsg] = mysql_dbexecute(dbc, updatestring);

		if ~isempty(emsg)
			disp(['Error Storing!! OH FUCK!!']); keyboard
		end		
	end

	times(i) = toc;
	if mod(i,floor(nempties/5))==0
		disp(['Approximately ' benrighttime(mean(times(1:i))*(nempties-i)) ...
			' remaining! -- We''ve so far written ' num2str(updatecounter) ' entries (/ ' num2str(i) ' possible)']);
	end
end

disp(['Well that seemed to go ok -- we wrote ' num2str(updatecounter) ' new entries / ' num2str(ntsgood*nmgood) ' to the database']);
if nempties > updatecounter % some were not written
	disp([num2str(nempties-updatecounter) ' entries were not written (recurring fatal errors) and remain awaiting calculation in the database']);
end


% 
% %% Check that outputs are ok
% disp(['Checking that we''re not messing things up in writing local outputs back to store...']);
% disp(['Re-retrieving TS_loc from database']);
% TS_loc_check = zeros(nts,nm);
% TS_loc_q_check = zeros(nts,nm);
% TS_loc_ct_check = zeros(nts,nm);
% times = zeros(nts,1);
% for i=1:nts
% 	selectstring = ['SELECT Output, Quality, CalculationTime FROM Results WHERE ts_id = ' num2str(ts_ids_keep(i)) ...
% 						' AND m_id IN (' m_ids_keep_string ' )'];
% 	[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% 	if ~isempty(emsg)
% 		disp(['Error retrieving outputs from STORE at ' num2str(ts_ids_keep(i)) ' -- Exiting']); keyboard
% 	end
% 	% Convert empty entries to NaNs
% 	empties = find(cellfun(@isempty,qrc));
% 	if ~isempty(empties), qrc(empties) = {NaN}; end
% 	
% 	TS_loc_check(i,:) = horzcat(qrc{:,1});
% 	TS_loc_q_check(i,:) = horzcat(qrc{:,2});
% 	TS_loc_ct_check(i,:) = horzcat(qrc{:,3});
% 	
% 	if mod(i,floor(nts/5))==0
% 		disp(['Approximately ' benrighttime(mean(times(1:i))*(nts-i)) ' remaining!']);
% 	end
% end
% 
% % Where there are values in TS_loc_check, we want to make sure that none of these are different.
% % in the new one -- this would be weird, because they shouldn't have been touched at all (could only happen
% % if another machine has calculated them in the interim...)
% isfisd = isfinite(TS_loc_check);
% % only check up to single precision
% notgoods = (single(TS_loc_check(isfisd)) ~= single(TS_loc(isfisd))); % but now storing at lower precision=, as singles, so "=" need not be exact
% 
% 
% if any(notgoods)
% 	
%     ncnf = (single(TS_loc_check)~=single(TS_loc) & isfinite(TS_loc_check)); % not consistent and finite in STORE
%     
%     [ai,bi] = find(ncnf);
%     uai = unique(ai); % rows that contain inconsistent entries
%     ubi = unique(bi); % columns that contain inconsistent entries
% 	Nincts = length(unique(ai)); % number of inconsistent time series
% 	Nincm = length(unique(bi)); % number of inconsistent operations
% 
% 	disp([num2str(sum(notgoods)) ' results are inconsistent']);
% 	disp(['There are ' num2str(Nincts) ' inconsistent time series']);
% 	disp(['There are ' num2str(Nincm) ' inconsistent metrics']);
% 	keyboard
% 	
% % 	% Likely a metric problem: remove all bad metrics
% % 	disp(['Trimming down to ' num2str(length(setxor(1:length(mgoodi),ubi))) ' metrics, if that''s ok?']);
% %     input('Control C now or forever hold your peace')
% %     
% % %     mgoodi=mgoodi(setxor(1:end,ubi));
% % %     keyboard
% % 	mgoodi(ubi)=[];
% else
% 	disp('None of the new results are inconsistent with existing results in the storage file. This is a good thing.');
% end
% 
% % % This next little section is probably not necessary; just to check I've done the right thing above...
% % STORE_data_tmp = STORE_data(tsmap(tsgoodi),mmap(mgoodi));
% % TS_loc_tmp = TS_loc(tsgoodi,mgoodi);
% % isfisd = isfinite(STORE_data_tmp);
% % if any(STORE_data_tmp(isfisd)~=TS_loc_tmp(isfisd))
% % 	beep, keyboard, return
% % end
% 
% %% Do the saving
% % disp(['Good part of TS_loc_AGG has ' num2str(sum(isfinite(TS_loc(:)))-sum(isfinite(STORE_data_tmp(:)))) ' excess real entries that will be written to store']);
% disp(['About to write the output of ' num2str(ntsgood) ' time series and ' num2str(nmgood) ' metrics ' ...
% 		'from TS_loc, TS_loc_q, and TS_loc_cts to the database...'])
% 
% %% How many new real entries?
% newentries = isnan(TS_loc_q_check); % none of TS_loc_q are NaNs -- there is something to be said for all of these now.
% % TS_loc_q_check is NaN but TS_loc is 0
% nnewrealentries = sum(sum(newentries & TS_loc_q==0));
% disp(['(*****) There are ' num2str(nnewrealentries) ' good, ' num2str(sum(sum(newentries))) ' new / ' num2str(ntsgood*nmgood) ' entries! Write away, my son! (******)']);
% 
% %% Do the writing
% times = zeros(nts,1);
% for i=1:nts
% 	tic
% 	for j=1:nm
% 		% if TS_loc(i,j)~=TS_loc_check(i,j) || TS_loc_q(i,j)~=TS_loc_q_check(i,j) || TS_loc_ct(i,j)~=TS_loc_ct_check(i,j);
% 		if newentries(i,j)
% 			% row has changed from calculation -- update with new values (updates LastModified automatically)
% 
% 			% Calcultion time -- make NULL if NaN
% 			if isnan(TS_loc_ct(i,j))
% 				string_ct = 'NULL';
% 			else
% 				string_ct = num2str(TS_loc_ct(i,j));
% 			end
% 
% 			updatestring = ['UPDATE Results SET ' ...
% 								'Output = ' num2str(TS_loc(i,j),'%19.17g') ', ' ...
% 								'Quality = ' num2str(TS_loc_q(i,j)) ', ' ...
% 								'CalculationTime = ' string_ct ', ' ...
% 								'LastModified = NOW() ' ...
% 								'WHERE ts_id = ' num2str(ts_ids_keep(i)) ' AND m_id = ' num2str(m_ids_keep(j))];
% 
% 			% 'CalculationTime = ' num2str(TS_loc_ct(i,j)) ', ' ...
% 		    [rs,emsg] = mysql_dbexecute(dbc, updatestring);
% 			if ~isempty(emsg)
% 				disp(['Error Storing!! OH FUCK!!']); keyboard
% 			end
% 		end
% 	end
% 	times(i) = toc;
% 	if mod(i,floor(nts/10))==0
% 		disp(['Approximately ' benrighttime(mean(times(1:i))*(nts-i)) ' remaining!']);
% 	end
% end

disp('All updates stored back in database! Hurrah!!');

SQL_closedatabase(dbc) % close database connection
disp(['Database ' char(dbc) ' closed']);

%% Write a log file
if tolog
	% disp('<TS_agglomerate>: Logging progess to file');
	% fn = ['TS_agglomerate_' datestr(now,30) '.log'];
	% fid = fopen(fn,'w','n');
	% disp(['Log file created: ' fn]);

	disp(['Logging to file...']);
	fn = ['TS_agglomerate_' datestr(now,30) '.log'];
	fid = fopen(fn,'w','n');
	disp(['Log file created: ' fn]);

	fprintf(fid, '%s\n', ['Updated ' num2str(length(tsgoodi)) ' time series: ']);
	for i=1:length(tsgoodi)
		fprintf(fid, '%s\n',tsf{tsgoodi(i)});
	end

	fprintf(fid, '\n\n\n\n\n%s\n', '******************************');
	fprintf(fid, '%s\n', ['Updated ' num2str(length(mgoodi)) ' metrics: ']);
	for i=1:length(mgoodi)
		fprintf(fid, '%s\n',mlab{mgoodi(i)});
	end
	fclose(fid);
end


end