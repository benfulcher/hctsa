% TSQ_prepared
% 
% This function retreives data from the mySQL database for subsequent analysis
% in Matlab. It takes as input a set of constraints on the time series and
% operations to include, and outputs the relevant subsection of the data matrix
% and associated metadata in HCTSA_loc
% 
%---INPUTS:
%--ts_ids: a vector of ts_ids to retrieve from the mySQL database.
%--op_ids: a vector of op_ids to retrieve from the mySQL database.
%--RetrieveWhatEntries: can be one of the following options:
%        (i) 'null': Retrieve null entries from the database (e.g., for
%                    calculating)
%        (ii) 'all': Retrieve all entries from the database (e.g., for
%                    analyzing)
%        (iii) 'error': Retrieve previous errors stored in the database (e.g.,
%                       for re-evaluating previous errors)
%--RetrieveWhatData: can be one of the following options:
%       (i) 'all': Retrieves Data, CalcTimes, and Quality labels
%       (ii) 'nocalctime': Retrieves calculated values and quality labels
%       (iii) 'outputs': Just the calculated values
%       (iv) 'quality': Just the quality labels
% 
%---OUTPUT:
%--DidWrite [opt] is 1 if Writes new file HCTSA_loc.mat
% 
% Other outputs are to the file HCTSA_loc.mat contains
%--TS_DataMat, contains the data
%--TS_Quality, contains the quality codes
%--TS_CalcTime, contains the calculation times [if RetrieveWhatData = 'all']
%--TimeSeries, contains metadata about the time series, and the time-series data
%--Operations, contains metadata about the operations
%--MasterOperations, contains metadata about the implicatedmaster operations
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

function DidWrite = TSQ_prepared(ts_ids, op_ids, RetrieveWhatEntries, RetrieveWhatData)
    
% Until it actually writes, set DidWrite = 0
DidWrite = 0;

%% Check inputs and set defaults
if nargin < 2
	error('You must provide at least two inputs!');
end

% RetrieveWhatEntries
if nargin < 3 || isempty(RetrieveWhatEntries)
	RetrieveWhatEntries = 'all'; % retrieve full sets of things, not just the empty entries in the database
    fprintf(1,'Retrieving ALL elements from the database in the specified ts_id and op_id ranges by default.\n')
end
RetrieveWhatEntriesCanBe = {'null','all','error'};
if ~ischar(RetrieveWhatEntries) || ~ismember(RetrieveWhatEntries,RetrieveWhatEntriesCanBe)
    error(['The third input to TSQ_prepared must specify what type of entries to retrieve from ' ...
                'the database, one of the following: %s'],BF_cat(RetrieveWhatEntriesCanBe))
end

% RetrieveWhatData
if nargin < 4 || isempty(RetrieveWhatData)
    RetrieveWhatData = 'nocalctime'; % Just get data and quality labels
end
RetrieveWhatDataCanBe = {'all','nocalctime','outputs','quality'};
if ~ischar(RetrieveWhatData) || ~ismember(RetrieveWhatData,RetrieveWhatDataCanBe)
    error(['The fourth input to TSQ_prepared must specify what data to retrieve from the database, ' ...
                'one of the following: %s'],BF_cat(RetrieveWhatEntriesCanBe))
end

%% Get going
% Retrieve from the database over many iterations, retrieving data from a given
% number of time series at each iteration
% Could do one big query and then reform to a matrix, but I'll do it row-by-row
% This is faster than doing one big query

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
if length(tsids_db) < nts % Actually there are fewer time series in the database than requested
    if (length(tsids_db) == 0) % Now there are no time series to retrieve
        fprintf(1,'None of the %u specified time series exist in ''%s''\n',nts,dbname)
        SQL_closedatabase(dbc); return % Close the database connection before returning
    end
    fprintf(1,['%u specified time series do not exist in ''%s'', retrieving' ...
                    ' the remaining %u\n'],nts-length(tsids_db),dbname,length(tsids_db))
    ts_ids = tsids_db;
    ts_ids_string = BF_cat(ts_ids,',');
    nts = length(ts_ids);
end
if length(opids_db) < nops % actually there are fewer operations in the database
    if (length(opids_db) == 0) % now there are no operations to retrieve
        fprintf(1,'None of the %u specified operations exist in ''%s''\n',nops,dbname)
        SQL_closedatabase(dbc); return % Close the database connection before returning
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
switch RetrieveWhatData
case 'all'
    TS_DataMat = ones(nts,nops)*Inf;  % Outputs
    TS_Quality = ones(nts,nops)*Inf;  % Quality labels
    TS_CalcTime = ones(nts,nops)*Inf; % Calculation times
case 'nocalctime'
    TS_DataMat = ones(nts,nops)*Inf;  % Outputs
    TS_Quality = ones(nts,nops)*Inf;  % Quality labels
case 'outputs'
    TS_DataMat = ones(nts,nops)*Inf;  % Outputs
case 'quality'
    TS_Quality = ones(nts,nops)*Inf;  % Quality labels
end

% Set BundleSize for retrieving
BundleSize = min(nts,1); % Retrieve data for this many time series per database query
                         % either 2 at a time, or if less, however many time series there are

%% Display information
switch RetrieveWhatEntries
case 'all' 
    fprintf(1,['Retrieving all elements from the database (in groups of %u time series ' ...
                'per database query). Please be patient...\n'],BundleSize);
case 'null'
    fprintf(1,['Retrieving NULL elements from the database (in groups of %u time series ' ...
                'per database query). Please be patient...\n'],BundleSize);
case 'error'
    fprintf(1,['Retrieving error elements from the database (in groups of %u time series ' ...
                'per database query). Please be patient...\n'],BundleSize);
end

Bundles = (1:BundleSize:length(ts_ids)); % Bundles of time series for each iteration
NumIterations = length(Bundles);         % Number of iterations of the loop
IterationTimes = zeros(NumIterations,1); % Record the time taken for each iteration
DidRetrieve = zeros(NumIterations,1);    % Keep track of whether data was retrieved at each iteration

% First select the data:
for i = 1:NumIterations
    
    IterationTimer = tic; % Time each iteration using IterationTimer

    ii = (Bundles(i):1:min(Bundles(i)+BundleSize-1,length(ts_ids))); % Indicies for this iteration
    ts_ids_now = ts_ids(ii); % Range of ts_ids retrieved in this iteration
    
    % Start piecing together the mySQL SELECT command:
    switch RetrieveWhatData
    case 'all'
        SelectWhat = 'SELECT ts_id, op_id, Output, QualityCode, CalculationTime FROM Results';
    case 'nocalctime' % Do not retrieve calculation time results
        SelectWhat = 'SELECT ts_id, op_id, Output, QualityCode FROM Results';
    case 'outputs'
        SelectWhat = 'SELECT ts_id, op_id, Output FROM Results';
    case 'quality'
        SelectWhat = 'SELECT ts_id, op_id, QualityCode FROM Results';
    end
    
    BaseString = sprintf('%s WHERE ts_id IN (%s) AND op_id IN (%s)',SelectWhat, ...
                                            BF_cat(ts_ids_now,','),op_ids_string);
    
    switch RetrieveWhatEntries
    case 'all'
        SelectString = BaseString;
    case 'null'
    	SelectString = sprintf('%s AND QualityCode IS NULL',BaseString);
    case 'error'
    	SelectString = sprintf('%s AND QualityCode = 1',BaseString);
    end
    
    % Do the retrieval
    % DatabaseTimer = tic;
	[qrc, ~, ~, emsg] = mysql_dbquery(dbc,SelectString); % Retrieve this bundle of data from the database
    % fprintf(1,'Database query for %u time series took %s\n',BundleSize,BF_thetime(toc(DatabaseTimer)));
    
    if ~isempty(emsg)
        error('Error retrieving outputs from %s!!! :(\n%s',dbname,emsg);
    end
    
    % Check results look ok:
    if (size(qrc) == 0) % There are no entries in Results that match the requested conditions
        fprintf(1,'No data to retrieve for ts_id = %s\n',BF_cat(ts_ids_now,','));
        % Leave local files TS_DataMat, TS_Quality (and TS_CalcTime) as Inf
        
    else
        % Entries need to be written to local matrices
        % Set DidRetrieve = 1 for this iteration
        DidRetrieve(i) = 1;
        
    	% Convert empty entries to NaNs
    	qrc(cellfun(@isempty,qrc)) = {NaN};
    
        % Put results from database into rows of local matrix
        if strcmp(RetrieveWhatEntries,'all') % easy in this case
            switch RetrieveWhatData
            case 'all'
                TS_DataMat(ii,:) = reshape(vertcat(qrc{:,3}),length(ii),nops);
                TS_Quality(ii,:) = reshape(vertcat(qrc{:,4}),length(ii),nops);
                TS_CalcTime(ii,:) = reshape(vertcat(qrc{:,5}),length(ii),nops);
            case 'nocalctime'
                TS_DataMat(ii,:) = reshape(vertcat(qrc{:,3}),length(ii),nops);
                TS_Quality(ii,:) = reshape(vertcat(qrc{:,4}),length(ii),nops);
            case 'outputs'
                TS_DataMat(ii,:) = reshape(vertcat(qrc{:,3}),length(ii),nops);
            case 'quality'
                TS_Quality(ii,:) = reshape(vertcat(qrc{:,4}),length(ii),nops);
            end
        else
            % We need to match retrieved indicies to the local indicies
            ix = arrayfun(@(x)find(ts_ids_now == x,1),vertcat(qrc{:,1}));
            iy = arrayfun(@(x)find(op_ids == x,1),vertcat(qrc{:,2}));
            switch RetrieveWhatData
            case 'all'
                for k = 1:length(ix) % fill it one entry at a time
                    TS_DataMat(ix(k),iy(k)) = qrc{k,3};
                    TS_Quality(ix(k),iy(k)) = qrc{k,4};
                    TS_CalcTime(ix(k),iy(k)) = qrc{k,5};
                end
            case 'nocalctime'
                for k = 1:length(ix) % fill it one entry at a time
                    TS_DataMat(ix(k),iy(k)) = qrc{k,3};
                    TS_Quality(ix(k),iy(k)) = qrc{k,4};
                end
            case 'outputs'
                for k = 1:length(ix) % fill it one entry at a time
                    TS_DataMat(ix(k),iy(k)) = qrc{k,3};
                end
            case 'quality'
                for k = 1:length(ix) % fill it one entry at a time
                    TS_Quality(ix(k),iy(k)) = qrc{k,4};
                end
            end
        end
    end
    
    % Note time taken for this iteration, and periodically display indication of time remaining
	IterationTimes(i) = toc(IterationTimer);
    if (i==1) % Give an initial indication of time after the first iteration
        fprintf(1,['Based on the first retrieval, this is taking ' ...
                'approximately %s per time series...\n'],BF_thetime(IterationTimes(1)/BundleSize));
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(IterationTimes(1)*(NumIterations-1)));
    elseif (mod(i,floor(NumIterations/10))==0) % Tell us the time remaining 10 times across the total retrieval
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(mean(IterationTimes(1:i))*(NumIterations-i)));
	end
end

% Finished retrieving from the database!
if any(DidRetrieve)
    fprintf(1,'Retrieved data from %s in %u iterations in %s.\n',dbname,NumIterations,BF_thetime(sum(IterationTimes)));
else
    fprintf(1,'Over %u iterations, no data was retrieved from %s.\nNot writing any data to file.\n',NumIterations,dbname);
    SQL_closedatabase(dbc); return
end
	
if ismember(RetrieveWhatEntries,{'null','error'})    
    % We only want to keep rows and columns with (NaNs for 'null' or errors for 'error') in them...
    switch RetrieveWhatEntries
    case 'null'
        keepme = isnan(TS_DataMat); % NULLs in database
        fprintf(1,['Filtering so that local files contain rows/columns containing at least ' ...
                                'one entry that was NULL in the database.\n']);
    case 'error'
        keepme = (TS_Quality == 1); % Error codes in database
        fprintf(1,['Filtering so that local files contain rows/columns containing at least ' ...
                                'one entry that was an error in the database.\n']);
    end
    
	% Time series
    keepi = (sum(keepme,2) > 0); % there is at least one entry to calculate in this row
    if sum(keepi) == 0
    	fprintf(1,'After filtering, there are no time series remaining! Exiting...\n');
        SQL_closedatabase(dbc); return % Close the database connection, then exit
	elseif sum(keepi) < nts
		fprintf(1,'Cutting down from %u to %u time series\n',nts,sum(keepi));
		ts_ids = ts_ids(keepi); nts = length(ts_ids);
		ts_ids_string = BF_cat(ts_ids,',');
        switch RetrieveWhatData
        case 'all'
    		TS_DataMat = TS_DataMat(keepi,:);
            TS_Quality = TS_Quality(keepi,:);
            TS_CalcTime = TS_CalcTime(keepi,:);
        case 'nocalctime'
    		TS_DataMat = TS_DataMat(keepi,:);
            TS_Quality = TS_Quality(keepi,:);
        case 'outputs'
            TS_DataMat = TS_DataMat(keepi,:);
        case 'quality'
            TS_Quality = TS_Quality(keepi,:);
        end
	end
	
	% Operations
    keepi = (sum(keepme,1) > 0); % there is at least one entry to calculate in this column
	if sum(keepi) == 0
    	fprintf(1,'After filtering, there are no operations remaining! Exiting...\n');
        SQL_closedatabase(dbc); return % Close the database connection, then exit
    elseif sum(keepi) < nops
		fprintf(1,'Cutting down from %u to %u operations\n',nops,sum(keepi));
		op_ids = op_ids(keepi); nops = length(op_ids);
		op_ids_string = BF_cat(op_ids,',');
        switch RetrieveWhatData
        case 'all'
    		TS_DataMat = TS_DataMat(:,keepi);
            TS_Quality = TS_Quality(:,keepi);
            TS_CalcTime = TS_CalcTime(:,keepi);
        case 'nocalctime' 
    		TS_DataMat = TS_DataMat(:,keepi);
            TS_Quality = TS_Quality(:,keepi);
        case 'outputs'
            TS_DataMat = TS_DataMat(:,keepi);
        case 'quality'
            TS_Quality = TS_Quality(:,keepi);
        end
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
save('HCTSA_loc.mat','TS_DataMat','TS_Quality','TimeSeries','Operations','MasterOperations','-v7.3');
switch RetrieveWhatData
case 'all'
    % Add outputs, quality labels, and calculation times
    save('HCTSA_loc.mat','TS_DataMat','TS_Quality','TS_CalcTime','-append')
case 'nocalctime'
    % Add outputs and quality labels
    save('HCTSA_loc.mat','TS_DataMat','TS_Quality','-append')
case 'outputs'
    % Add outputs
    save('HCTSA_loc.mat','TS_DataMat','-append')
case 'quality'
    % Add quality labels
    save('HCTSA_loc.mat','TS_Quality','-append')
end

fprintf(1,' Done.\n');
DidWrite = 1;

% Display how many entries need to be calculated
tocalculate = sum(isnan(TS_Quality(:)) | TS_Quality(:)==1);
fprintf(1,'There are %u entries (=%4.2f%%) to calculate in the data matrix (%ux%u).\n', ...
                                    tocalculate,tocalculate/nops/nts*100,nts,nops);

end