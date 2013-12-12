% TSQ_prepared
% 
% This function retreives data from the mySQL database for subsequent analysis
% in Matlab. It takes as input a set of constraints on the time series and
% operations to include, and outputs the relevant subsection of the data matrix
% and associated metadata in HCTSA_loc
% 
% ---OUTPUT to the file HCTSA_loc.mat contains
% (*) TS_DataMat, contains the data
% (*) TS_CalcTime, contains the calculation times
% (*) TS_Quality, contains the quality codes
% (*) TimeSeries, contains metadata about the time series, and the time-series data
% (*) Operations, contains metadata about the operations
% (*) MasterOperations, contains metadata about the implicatedmaster operations
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TSQ_prepared(ts_ids, op_ids, RetrieveWhat, brawninputs, doinone)
%% Check inputs -- set defaults
if nargin < 2
	error('You must provide at least two inputs!');
end
if nargin < 3
	RetrieveWhat = 'all'; % retrieve full sets of things, not just the empty entries in the database
    fprintf(1,'Retrieving ALL elements from the database in the specified ranges by default.\n')
end
RetrieveWhatcanbe = {'null','all','error'};
if ~ischar(RetrieveWhat) || ~ismember(RetrieveWhat,RetrieveWhatcanbe)
    error('The rhird input to TSQ_prepared must specify what to retrieve, one of the following: %s',BF_cat(RetrieveWhatcanbe))
end
if nargin < 4 || isempty(brawninputs)
	brawninputs = [1, 1]; % log and parallelize -- only relevant if RetrieveWhat isn't empty
end
if nargin < 5 || isempty(doinone)
	doinone = 0; % make seperate connections so as not to overwhelm java heap space
end

%% METHOD 1: Get entries of results table for local Matlab matrix (1 row/column at a time)
% we have ts_ids and op_ids
% To put in a matrix with rows (time series) and columns (metrics)
% Could do one big query and then reform to a matrix, but I'll do it row-by-row
% In fact this is faster for some reason than doing a big query (method 2)

% Make sure ts_ids and op_ids are column vectors
if size(ts_ids,2) > size(ts_ids,1), ts_ids = ts_ids'; end
if size(op_ids,2) > size(op_ids,1), op_ids = op_ids'; end
% Sort ids ascending
ts_ids = sort(ts_ids,'ascend');
op_ids = sort(op_ids,'ascend');
% Write a comma-delimited string of ids
ts_ids_string = BF_cat(ts_ids,',');
op_ids_string = BF_cat(op_ids,',');
% Count the number of time series and operations
nts = length(ts_ids); nops = length(op_ids);

if (nts == 0)
	error('Oops. There''s nothing to do! No time series to retrieve!\n');
elseif (nops == 0)
	error('Oops. There''s nothing to do! No operations to retrieve!\n');
end

% Open database connection
[dbc, dbname] = SQL_opendatabase;

% First refine the set of time series and operations to those that actually exist in the database
opids_db = mysql_dbquery(dbc,sprintf('SELECT op_id FROM Operations WHERE op_id IN (%s)',op_ids_string));
opids_db = vertcat(opids_db{:});
tsids_db = mysql_dbquery(dbc,sprintf('SELECT ts_id FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string));
tsids_db = vertcat(tsids_db{:});
if length(tsids_db) < nts % actually there are fewer time series in the database
    if (length(tsids_db) == 0) % now there are no time series to retrieve
        error('None of the %u specified time series exist in ''%s''',nts-length(tsids_db),dbname)
    end
    fprintf(1,['%u specified time series do not exist in ''%s'', retrieving' ...
                    ' the remaining %u\n'],nts-length(tsids_db),dbname,length(tsids_db))
    ts_ids = tsids_db;
    ts_ids_string = BF_cat(ts_ids,',');
    nts = length(ts_ids);
end
if length(opids_db) < nops % actually there are fewer operations in the database
    if (length(opids_db) == 0) % now there are no operations to retrieve
        error('None of the %u specified operations exist in ''%s''',nops-length(opids_db),dbname)
    end
    fprintf(1,['%u specified operations do not exist in ''%s'', retrieving' ...
                    ' the remaining %u\n'],nops-length(opids_db),dbname,length(opids_db))
    op_ids = opids_db;
    op_ids_string = BF_cat(op_ids,',');
    nops = length(op_ids);
end

% Tell me about it
fprintf(1,'We have %u time series and %u operations to retrieve from %s.\n',nts,nops,dbname);
fprintf(1,'Filling and saving to local Matlab file HCTSA_loc.mat from the Results table of %s.\n',dbname);

%% Intialize matrices

% Start as Infs to distinguish unwritten entries after database pull
TS_DataMat = ones(nts,nops)*Inf; % outputs
TS_CalcTime = ones(nts,nops)*Inf; % calculation times
TS_Quality = ones(nts,nops)*Inf; % output quality label

% Set bundlesize for retrieving
bundlesize = min(nts,5); % retrieve information about this many time series per database query:
                         % either 5 at a time, or if less, however many time series there are

%% Provide user information
switch RetrieveWhat
case 'all' 
    fprintf(1,'Retrieving all elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
case 'null'
    fprintf(1,'Retrieving NULL elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
case 'error'
    fprintf(1,'Retrieving error elements from the database (in groups of %u time series, FYI). Please be patient...\n',bundlesize);
end

bundles = (1:bundlesize:length(ts_ids));
nits = length(bundles); % number of iterations of the loop
times = zeros(nits,1); % record time for each iteration

% First select the data:
for i = 1:nits
    ittic = tic; % start timing the iteration

    ii = (bundles(i):1:min(bundles(i)+bundlesize-1,length(ts_ids))); % indicies for this iteration
    ts_ids_now = ts_ids(ii); % a range of ts_ids to retrieve in this iteration
    basestring = sprintf(['SELECT ts_id, op_id, Output, CalculationTime, QualityCode FROM Results WHERE ' ...
                    	'ts_id IN (%s) AND op_id IN (%s)'],BF_cat(ts_ids_now,','),op_ids_string);
    switch RetrieveWhat
    case 'all'
        SelectString = basestring;
    case 'null'
    	SelectString = [basestring, ' AND QualityCode IS NULL'];
    case 'error'
    	SelectString = [basestring, ' AND QualityCode = 1'];
    end
    
	[qrc, ~, ~, emsg] = mysql_dbquery(dbc,SelectString); % retrieve the bundlesize from the database
    
    % Check results look ok:
    if ~isempty(emsg)
        fprintf(1,'Error retrieving outputs from %s???\n',dbname);
        fprintf(1,'%s\n',emsg)
        keyboard
    end
    
    if (size(qrc) == 0)
        fprintf(1,'No data to retrieve for ts_id = (%s)\n',BF_cat(ts_ids_now,','));
    end
    
	% Convert empty entries to NaNs
	qrc(cellfun(@isempty,qrc)) = {NaN};
    
    % Fill rows of local matrix
    if strcmp(RetrieveWhat,'all') % easy in this case
        TS_DataMat(ii,:) = reshape(vertcat(qrc{:,3}),length(ii),nops);
        TS_CalcTime(ii,:) = reshape(vertcat(qrc{:,4}),length(ii),nops);
        TS_Quality(ii,:) = reshape(vertcat(qrc{:,5}),length(ii),nops);
    else
        % We need to match retrieved indicies to the local indicies
        ix = arrayfun(@(x)find(ts_ids_now == x,1),vertcat(qrc{:,1}));
        iy = arrayfun(@(x)find(op_ids == x,1),vertcat(qrc{:,2}));
        for k = 1:length(ix) % fill it one entry at a time
            TS_DataMat(ix(k),iy(k)) = qrc{k,3};
            TS_CalcTime(ix(k),iy(k)) = qrc{k,4};
            TS_Quality(ix(k),iy(k)) = qrc{k,4};
        end
    end
    
    % Note time taken for this iteration
	times(i) = toc(ittic);
	if mod(i,floor(nits/10)) == 0 % tell the user 10 times
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(mean(times(1:i))*(nts-i)));
	end
end
fprintf(1,'Local files filled from %s in %u iteration(s). Took %s altogether.\n',dbname,nits,BF_thetime(sum(times)));
	
if ismember(RetrieveWhat,{'null','error'})    
    % We only want to keep rows and columns with (e.g., NaNs) in them...
    switch RetrieveWhat
    case 'null'
        keepme = isnan(TS_DataMat); % NULLs in database
    	fprintf(1,'Filtering so that local files contain rows/columns containing at least one NULL entry\n');
    case 'error'
        keepme = (TS_Quality == 1); % Error codes in database
    	fprintf(1,'Filtering so that local files contain rows/columns containing at least one error entry\n');
    end
    
	% Time series
    keepi = (sum(keepme,2) > 0); % there is at least one entry to calculate in this row
    if sum(keepi) == 0
    	fprintf(1,'After filtering, there are no time series remaining! Exiting...\n'); return
	elseif sum(keepi) < nts
		fprintf(1,'Cutting down from %u to %u time series\n',nts,sum(keepi));
		ts_ids = ts_ids(keepi); nts = length(ts_ids);
		ts_ids_string = BF_cat(ts_ids,',');
		TS_DataMat = TS_DataMat(keepi,:); TS_CalcTime = TS_CalcTime(keepi,:); TS_Quality = TS_Quality(keepi,:); % update local data stores
	end
	
	% Operations
    keepi = (sum(keepme,1) > 0); % there is at least one entry to calculate in this column
	if sum(keepi) == 0
    	fprintf(1,'After filtering, there are no operations remaining! Exiting...\n'); return
    elseif sum(keepi) < nops
		fprintf(1,'Cutting down from %u to %u operations\n',nops,sum(keepi));
		op_ids = op_ids(keepi); nops = length(op_ids);
		op_ids_string = BF_cat(op_ids,',');
		TS_DataMat = TS_DataMat(:,keepi); TS_CalcTime = TS_CalcTime(:,keepi); TS_Quality = TS_Quality(:,keepi);
	end    
end

%% Fill Metadata

% 1. Retrieve Time Series Metadata
SelectString = sprintf('SELECT FileName, Keywords, Length, Data FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string);
[tsinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
% Convert to a structure array, TimeSeries, containing metadata for all time series
tsinfo = [num2cell(ts_ids),tsinfo];
% Define inline functions to convert time-series data text to a vector of floats:
ScanCommas = @(x) textscan(x,'%f','Delimiter',',');
TakeFirstCell = @(x) x{1};
tsinfo(:,end) = cellfun(@(x) TakeFirstCell(ScanCommas(x)),tsinfo(:,end),'UniformOutput',0); % Do the conversion
TimeSeries = cell2struct(tsinfo',{'ID','FileName','Keywords','Length','Data'}); % Convert to structure array

% 2. Retrieve Operation Metadata
SelectString = sprintf('SELECT OpName, Keywords, Code, mop_id FROM Operations WHERE op_id IN (%s)',op_ids_string);
[opinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
opinfo = [num2cell(op_ids), opinfo]; % add op_ids
Operations = cell2struct(opinfo',{'ID','Name','Keywords','CodeString','MasterID'});

% 3. Retrieve Master Operation Metadata
% (i) Which masters are implicated?
SelectString = ['SELECT mop_id, MasterLabel, MasterCode FROM MasterOperations WHERE mop_id IN ' ...
        				'(' BF_cat(unique([Operations.MasterID]),',') ')'];
[masterinfo,~,~,emsg] = mysql_dbquery(dbc,SelectString);
if ~isempty(emsg)
    fprintf(1,'Error retrieving Master information...\n'); keyboard
else
    MasterOperations = cell2struct(masterinfo',{'ID','Label','Code'});
end

% Close database connection
SQL_closedatabase(dbc)

% Save to HCTSA_loc.mat
fprintf(1,'Saving local versions of the data...');
save('HCTSA_loc.mat','TS_DataMat','TS_CalcTime','TS_Quality','TimeSeries','Operations','MasterOperations','-v7.3');
fprintf(1,' Done.\n');

% Display how many entries need to be calculated
tocalculate = sum(isnan(TS_Quality(:)) | TS_Quality(:)==1);
fprintf(1,'There are %u entries (=%4.2f%%) to calculate in the data matrix (%ux%u).\n', ...
                                    tocalculate,tocalculate/nops/nts*100,nts,nops);

end