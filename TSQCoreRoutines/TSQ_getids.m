function ids = TSQ_getids(morts,lenr,kyesc,kno,pcalcr,idr,howtolimit,dbname)

% This function takes as input a set of constraints on the time series and metrics to include
% and output the relevant ts_ids / m_ids

% HISTORY:
% 11/1/2010: Ben Fulcher ~ Imported from TSQ_prepared and rehauled for this purpose
% 10/4/2010: Ben Fulcher -- added howtolimit as new input

%% INPUTS:
% 1) morts: specifies 'mets' for metrics or 'ts' for time series
% 2) lenr (ts): contrain included time series by length [minimum_length maximum_length] (1x2 vector)
% 2) lenr (mets): takes roll of masterpull: whether (1) or not (0) to pull
%                   in other master-linked operations
% 3) kyesc: constrain included time series by keyword -- nx2 cell
% 			(i) keyword
% 			(ii) how many of that keyword (0==all available)
% 4) kno: keywords NOT to include (cell of strings)
% 5) pcalcr: range of percentagecalc to consider (empty = all, default)
% 6) idr: range of possible ids to restrict the search to (empty = don't
%                       restrict -- default)
% 7) howtolimit: can select one of:
%                   'rand' (limit cases by random)
%                   'pcmax' (limit cases by maximum percentage calculated)
%                   'pcmin' (limit cases by minimum percentage calculated)
% 8) dbc: can specify database; otherwise opens the default specified in
%         SQL_opendatabase

%% OUTPUTS
% ids: a vector or either ts_ids or m_ids that match the input constraints


%%% FOREPLAY
%% Check inputs -- set defaults
% if nargin < 5; disp('You must provide 5 inputs! And no, I''m not asking nicely!!'); return; end

if strcmp(morts,'ts')
	if isempty(lenr)
		lenr = [200, 20000];
		fprintf(1,'Setting default length constraints: %u--%u\n',lenr(1),lenr(2))
	end
elseif strcmp(morts,'mets')
	masterpull = lenr;
	if isempty(masterpull)
		masterpull = 1; % retrieves other outputs of all master functions implicated in the range
		fprintf(1,'Pulling in all pointers by default');
	end
else
	disp('First input must be either ''ts'' or ''mets''. I''m outraged. Exiting.');
    return
end

if nargin < 4
    kno = {};
end

if nargin < 5
	pcalcr = [];
end
if nargin < 6
	idr = [];
end
if nargin < 7 || isempty(howtolimit)
    howtolimit = 'pcmax';
end
if nargin < 8
	dbname = '';
end

% seperate into keywords to include: kyes
% and the number of that keyword to include: kyesn
if ~isempty(kyesc)
    kyes = kyesc(:,1); kyesn = vertcat(kyesc{:,2});
else
    kyes = {}; kyesn = [];
end


% if tolog
% 	fn = ['TSQ_prepared_' datestr(now,30) '.log'];
% 	flog = fopen(fn,'w','n');
% 	disp(['Logging progress to ' fn]);
% 	fprintf(flog,'%s\n',['Welcome! This is your friendly TSQ_prepared log']);
% 	fprintf(flog,'%s\n',['Subsetting and checking performed at ' datestr(now)]);
% end

%% Open connection to mySQL database
dbc = SQL_opendatabase(dbname);

%%% TIME SERIES
if strcmp(morts,'ts')

	%% Filter time series by keyword
	% keywords as a cell: kyes; and how many of each to include: kyesn
	s = {'','','',''};
	
	% Extra qualifier to only look in a certain range
	% s{1} -- IDR
	if ~isempty(idr)
		s{1} = sprintf('ts_id BETWEEN %u AND %u',idr(1),idr(2));
	end

	% Extra qualifier to not include certain keywords -- this is performed seperately in the query for each keyword *to* include
	% s{2} -- kno
	if ~isempty(kno)
		s{2} = ['ts_id NOT IN (SELECT ts_id FROM TsKeywordsRelate WHERE tskw_id IN ' ...
							'(SELECT tskw_id FROM TimeSeriesKeywords WHERE Keyword IN (' bencat(kno,',','''') '))) '];
	end
	
	% Extra qualifier pcalcr to have PercentageCalculated only in a certain range
	% s{3} -- pcalcr
	if ~isempty(pcalcr)
		s{3} = sprintf('PercentageCalculated BETWEEN %f AND %f',pcalcr(1),pcalcr(2));
	end
	
	% length constraint in words
	s{4} = sprintf('TimeSeries.Length BETWEEN %u AND %u',lenr(1));
	
	% Combine these results into part of a mySQL query string
	conditions = '';
	for j = 1:length(s)
		if ~isempty(s{j})
			conditions = [conditions ' AND ' s{j}];
		end
	end

	%% kyes, keywords to include
	if ~isempty(kyes)
		C_tskwyes = cell(length(kyes),1);
		for i = 1:length(kyes)
			ncut = kyesn(i);

			% Now do the rest of the query at once: do the keyword matches and length constraints		
			if ncut == 0 % Get all keyword matches
				SelectString = ['SELECT ts_id FROM TimeSeries WHERE ' ...
								'ts_id IN (SELECT ts_id FROM TsKeywordsRelate WHERE tskw_id = '  ...
								'(SELECT tskw_id FROM TimeSeriesKeywords WHERE Keyword = ''' kyes{i} '''))' ...
								conditions];
			else % constrain to some number (using the LIMIT command in mySQL)
                switch howtolimit
                    case 'pcmax'
                        limitme = [' ORDER BY PercentageCalculated DESC LIMIT ' num2str(ncut)];
                    case 'pcmin'
                        limitme = [' ORDER BY PercentageCalculated LIMIT ' num2str(ncut)];
                    case 'rand'
                        limitme = [' ORDER BY RAND() LIMIT ' num2str(ncut)];
                        % I know this is a slow implementation -- probably better
                        % to create the random numbers here to constrain
                end
                
				if isempty(kyes{i})
					SelectString = ['SELECT ts_id FROM TimeSeries WHERE ' ...
									conditions(6:end) limitme];

                else		
					% constrain by PercentageCalculated
					SelectString = ['SELECT ts_id FROM TimeSeries WHERE ' ...
									'ts_id IN (SELECT ts_id FROM TsKeywordsRelate WHERE tskw_id = '  ...
									'(SELECT tskw_id FROM TimeSeriesKeywords WHERE Keyword = ''' kyes{i} ''')) ' ...
									conditions limitme];
				end
			end

			% execute the select statement
			[ts_ids,~,~,emsg] = mysql_dbquery(dbc,SelectString);
			if isempty(ts_ids) && isempty(emsg)
				disp(['No time series matching constraints for ' kyes{i}]);
			elseif ~isempty(emsg)
				disp(['Error retrieving time series for ' kyes{i}]);
				disp(emsg)
				keyboard
			else
				C_tskwyes{i} = vertcat(ts_ids{:});
				ngot = length(C_tskwyes{i});
			
				if isempty(kyes{i})
					fprintf(1,'Found %u time series selected by %s\n',ngot,howtolimit);
				else
					fprintf(1,'Found %u time series with keyword ''%s''\n',ngot,kyes{i});
				end
			end
		end
	
		% Finally, now keep all that are left
		ts_ids_keep = vertcat(C_tskwyes{:});
		lbefore = length(ts_ids_keep);
		ts_ids_keep = unique(ts_ids_keep);
		lafter = length(ts_ids_keep);
		if lafter<lbefore
			fprintf(1,'We''ve lost %u time series to keyword overlap...\n',lbefore-lafter);
		end
	
	else % just use the other constraints
		% if isempty(kno)
			SelectString = ['SELECT ts_id FROM TimeSeries WHERE ' conditions(6:end)];
		% else
		% 	SelectString = ['SELECT ts_id FROM TimeSeries WHERE ' ss_tsl ' ' knoextrastring];
		% end
		[ts_ids,~,~,emsg] = mysql_dbquery(dbc,SelectString);
		if ~isempty(emsg)
			fprintf(1,'Database call failed\n'); disp(emsg), keyboard
		else
			ts_ids_keep = unique(vertcat(ts_ids{:}));
		end
	end

	% We now have ts_ids_keep -- the keyword and length-constrained time series
	nts = length(ts_ids_keep); % number of time series
	if nts == 0
		fprintf(1,'No Time Series found.\n');
		ids = [];
	else
		fprintf(1,'Time Series Filtered: %u\n',nts);
		ids = ts_ids_keep;
	end
	
else
	%%% Operations

	s = {'','',''};
	
	% Extra qualifier to only look in a certain range
	% s{1} -- IDR
	if ~isempty(idr)
		s{1} = ['m_id BETWEEN ' num2str(idr(1)) ' AND ' num2str(idr(2))];
	end

	% Extra qualifier to not include certain keywords -- this is performed seperately in the query for each keyword *to* include
	% s{2} -- kno
	if ~isempty(kno)
		s{2} = ['m_id NOT IN (SELECT m_id FROM OpKeywordsRelate WHERE mkw_id IN ' ...
							'(SELECT mkw_id FROM OperationKeywords WHERE Keyword IN (' bencat(kno,',','''') '))) '];
	end
	
	% Extra qualifier pcalcr to have PercentageCalculated only in a certain range
	% s{3} -- pcalcr
	if ~isempty(pcalcr)
		s{3} = ['PercentageCalculated BETWEEN ' num2str(pcalcr(1)) ' AND ' num2str(pcalcr(2))];
	end
	
	% Combine these results into part of a mySQL query string
	conditions='';
	for j = 1:length(s)
		if ~isempty(s{j})
			conditions = [conditions ' AND ' s{j}];
		end
	end
	
	% if ~isempty(kno)
	% 	mnoextrastring = ['AND m_id NOT IN (SELECT m_id FROM OpKeywordsRelate WHERE mkw_id IN ' ...
	% 						'(SELECT mkw_id FROM OperationKeywords WHERE Keyword IN (' bencat(mno,',','''') '))) '];
	% else
	% 	mnoextrastring = '';
	% end

	if ~isempty(kyes)
		C_kyes = cell(length(kyes),1);
		for i=1:length(kyes)
			ncut = kyesn(i);

			% Now do the rest of the query at once: do the keyword matches and length constraints
		
			if ncut==0 % Get all matches			
				SelectString = ['SELECT m_id FROM Operations WHERE ' ...
								'm_id IN (SELECT m_id FROM OpKeywordsRelate WHERE mkw_id = '  ...
								'(SELECT mkw_id FROM OperationKeywords WHERE Keyword = ''' kyes{i} '''))' ...
								conditions];
			else % constrain to some number (using the LIMIT command in mySQL)
                switch howtolimit
                    case 'pcmax'
                        limitme = [' ORDER BY PercentageCalculated DESC LIMIT ' num2str(ncut)];
                    case 'pcmin'
                        limitme = [' ORDER BY PercentageCalculated LIMIT ' num2str(ncut)];
                    case 'rand'
                        limitme = [' ORDER BY RAND() LIMIT ' num2str(ncut)];
                        % I know this is a slow implementation -- probably better
                        % to create the random numbers here to constrain
                end
                
                
				if isempty(kyes{i}) % empty keyword -- just constrain by number
					if isempty(conditions)
						SelectString = ['SELECT m_id FROM Operations' limitme];
					else
						SelectString = ['SELECT m_id FROM Operations WHERE ' conditions(6:end) ...
											limitme];
					end
				else
					SelectString = ['SELECT m_id FROM Operations WHERE ' ...
									'm_id IN (SELECT m_id FROM OpKeywordsRelate WHERE mkw_id = '  ...
									 '(SELECT mkw_id FROM OperationKeywords WHERE Keyword = ''' kyes{i} '''))' ...
									   conditions limitme];
				end
			end
		
			% Execute the query
			[m_ids,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
		
			if ~isempty(emsg)
				fprintf(1,'Error finding %s\n',kyes{i}); disp(emsg); keyboard
			end
			C_kyes{i} = vertcat(m_ids{:});
		
			ngot = length(C_kyes{i});
			if isempty(kyes{i})
				fprintf(1,'Found %u operations chosen by %s\n',ngot,howtolimit);
			else
				fprintf(1,'Found %u operations with keyword ''%s''\n',ngot,kyes{i});
			end
		end
	
		% Agglomerate all the bits
		m_ids_keep = vertcat(C_kyes{:});
		lbefore = length(m_ids_keep);
		m_ids_keep = unique(m_ids_keep);
		lafter = length(m_ids_keep);
		if lafter < lbefore
			fprintf(1,'We lost %u  to overlapping keywords!\n',lbefore-lafter);
		end
	
	else % just use the length/kno constraint
		if isempty(conditions)
			SelectString = 'SELECT m_id FROM Operations'; % include all operations -- exclude nothing
		else
			SelectString = ['SELECT m_id FROM Operations WHERE ' conditions(6:end)]; % (remove "AND ")
		end
		[m_ids,~,~,emsg] = mysql_dbquery(dbc,SelectString);
		if ~isempty(emsg)
			fprintf(1,'Database call failed\n'); disp(emsg); keyboard
		else
			m_ids_keep = unique(vertcat(m_ids{:}));
		end
	end


	%% Get other metrics which point to master functions which will be called anyway
	if masterpull && ~isempty(m_ids_keep)
		% Find implicated Master functions
		SelectString = ['SELECT m_id FROM MasterPointerRelate WHERE mop_id IN (SELECT DISTINCT mop_id FROM MasterPointerRelate WHERE m_id IN (' bencat(m_ids_keep,',') '))'];
		[newmids,~,~,emsg] = mysql_dbquery(dbc,SelectString);
		if isempty(emsg)
			if ~isempty(newmids) % there are some master functions implicated
				newmids = vertcat(newmids{:});
				nm = length(m_ids_keep); % number of metrics
				m_ids_keep = union(m_ids_keep,newmids); % include the new ones
				if length(m_ids_keep) > nm
					fprintf(1,'%u additional operations were included as implicated by existing master functions\n',length(m_ids_keep)-nm);
				end
			end
		else % an error
			fprintf(1,'Error retrieving the operations from implicated master functions\n'); disp(emsg); keyboard
		end
	end

	nm = length(m_ids_keep); % number of metrics
	if nm == 0
		fprintf(1,'No matching operations found.\n');
		ids = [];
	else
		fprintf(1,'Operations filtered: %u\n',nm);
		ids = m_ids_keep;
	end
end

% if tolog, fclose(flog), end
SQL_closedatabase(dbc)

end