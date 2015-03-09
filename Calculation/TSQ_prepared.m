% --------------------------------------------------------------------------
% TSQ_prepared
% --------------------------------------------------------------------------
% 
% This function retreives data from the mySQL database for subsequent analysis
% in Matlab. It takes as input a set of constraints on the time series and
% operations to include, and outputs the relevant subsection of the data matrix
% and associated metadata to HCTSA_loc.mat
% 
%---INPUTS:
%--ts_ids: a vector of ts_ids to retrieve from the mySQL database.
%--op_ids: a vector of op_ids to retrieve from the mySQL database.
%--retrieveWhatEntries: can be one of the following options:
%        (i) 'null': Retrieve null entries from the database (e.g., for
%                    calculating)
%        (ii) 'all': Retrieve all entries from the database (e.g., for
%                    analyzing)
%        (iii) 'error': Retrieve previous errors stored in the database (e.g.,
%                       for re-evaluating previous errors)
%--retrieveWhatData: can be one of the following options:
%       (i) 'all': Retrieves Data, CalcTimes, and Quality labels
%       (ii) 'nocalctime': Retrieves calculated values and quality labels
%       (iii) 'outputs': Just the calculated values
%       (iv) 'quality': Just the quality labels
% 
%---OUTPUT:
%--didWrite [opt] is 1 if Writes new file HCTSA_loc.mat
% 
% Other outputs are to the file HCTSA_loc.mat, which contains
%--TS_DataMat, contains the data
%--TS_Quality, contains the quality codes
%--TS_CalcTime, contains the calculation times [if retrieveWhatData = 'all']
%--TimeSeries, contains metadata about the time series, and the time-series data
%--Operations, contains metadata about the operations
%--MasterOperations, contains metadata about the implicatedmaster operations
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

function didWrite = TSQ_prepared(ts_ids, op_ids, retrieveWhatEntries, retrieveWhatData)
    
% Until it actually writes something, set the function output, didWrite = 0
didWrite = 0;

% --------------------------------------------------------------------------
%% Check inputs and set defaults
% --------------------------------------------------------------------------
if nargin < 2
	error('You must provide at least two inputs!');
end

% retrieveWhatEntries
if nargin < 3 || isempty(retrieveWhatEntries)
	retrieveWhatEntries = 'all'; % retrieve full sets of things, not just the empty entries in the database
    fprintf(1,'Retrieving ALL elements from the database in the specified ts_id and op_id ranges by default.\n')
end
retrieveWhatEntriesCanBe = {'null','all','error'};
if ~ischar(retrieveWhatEntries) || ~ismember(retrieveWhatEntries,retrieveWhatEntriesCanBe)
    error(['The third input to TSQ_prepared must specify what type of entries to retrieve from ' ...
                'the database, one of the following: %s'],BF_cat(retrieveWhatEntriesCanBe))
end

% retrieveWhatData
if nargin < 4 || isempty(retrieveWhatData)
    retrieveWhatData = 'nocalctime'; % Just get data and quality labels
end
retrieveWhatDataCanBe = {'all','nocalctime','outputs','quality'};
if ~ischar(retrieveWhatData) || ~ismember(retrieveWhatData,retrieveWhatDataCanBe)
    error(['The fourth input to TSQ_prepared must specify what data to retrieve from the database, ' ...
                'one of the following: %s'],BF_cat(retrieveWhatEntriesCanBe))
end

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------

% Make sure ts_ids and op_ids are column vectors:
if size(ts_ids,2) > size(ts_ids,1), ts_ids = ts_ids'; end
if size(op_ids,2) > size(op_ids,1), op_ids = op_ids'; end

% Sort ids ascending:
ts_ids = sort(ts_ids,'ascend');
op_ids = sort(op_ids,'ascend');

% Write a comma-delimited string of ids:
ts_ids_string = BF_cat(ts_ids,',');
op_ids_string = BF_cat(op_ids,',');

% Count the number of time series and operations
numTS = length(ts_ids);
numOps = length(op_ids);

if (numTS == 0)
	error('Oops. There''s nothing to do! No time series to retrieve!');
elseif (numOps == 0)
	error('Oops. There''s nothing to do! No operations to retrieve!');
end

% Open database connection
[dbc, dbname] = SQL_opendatabase;

% ------------------------------------------------------------------------------
% Refine the set of time series and operations to those that actually exist in the database
% ------------------------------------------------------------------------------
opids_db = mysql_dbquery(dbc,sprintf('SELECT op_id FROM Operations WHERE op_id IN (%s)',op_ids_string));
opids_db = vertcat(opids_db{:});
tsids_db = mysql_dbquery(dbc,sprintf('SELECT ts_id FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string));
tsids_db = vertcat(tsids_db{:});
if length(tsids_db) < numTS % Actually there are fewer time series in the database than requested
    if (length(tsids_db) == 0) % Now there are no time series to retrieve
        fprintf(1,'None of the %u specified time series exist in ''%s''\n',numTS,dbname)
        SQL_closedatabase(dbc); return % Close the database connection before returning
    end
    fprintf(1,['%u specified time series do not exist in ''%s'', retrieving' ...
                    ' the remaining %u\n'],numTS-length(tsids_db),dbname,length(tsids_db))
    ts_ids = tsids_db; % Will always be sorted in ascending order
    ts_ids_string = BF_cat(ts_ids,',');
    numTS = length(ts_ids);
end

if length(opids_db) < numOps % actually there are fewer operations in the database
    if (length(opids_db) == 0) % now there are no operations to retrieve
        fprintf(1,'None of the %u specified operations exist in ''%s''\n',numOps,dbname)
        SQL_closedatabase(dbc); return % Close the database connection before returning
    end
    fprintf(1,['%u specified operations do not exist in ''%s'', retrieving' ...
                    ' the remaining %u\n'],numOps-length(opids_db),dbname,length(opids_db))
    op_ids = opids_db; % Will always be sorted in ascending order
    op_ids_string = BF_cat(op_ids,',');
    numOps = length(op_ids);
end

% Tell me about it
fprintf(1,'We have %u time series and %u operations to retrieve from %s.\n',numTS,numOps,dbname);
fprintf(1,['Filling and saving to local Matlab file HCTSA_loc.mat from ' ...
                                'the Results table of %s.\n'],dbname);

% --------------------------------------------------------------------------
%% Intialize matrices
% --------------------------------------------------------------------------

% Initialize as Infs to distinguish unwritten entries after database retrieval
switch retrieveWhatData
case 'all'
    TS_DataMat = ones(numTS,numOps)*Inf;  % Outputs
    TS_Quality = ones(numTS,numOps)*Inf;  % Quality labels
    TS_CalcTime = ones(numTS,numOps)*Inf; % Calculation times
case 'nocalctime'
    TS_DataMat = ones(numTS,numOps)*Inf;  % Outputs
    TS_Quality = ones(numTS,numOps)*Inf;  % Quality labels
case 'outputs'
    TS_DataMat = ones(numTS,numOps)*Inf;  % Outputs
case 'quality'
    TS_Quality = ones(numTS,numOps)*Inf;  % Quality labels
end

% Display information to user:
switch retrieveWhatEntries
case 'all' 
    fprintf(1,['Retrieving all elements from the database (one time series ' ...
                'per database query). Please be patient...\n']);
case 'null'
    fprintf(1,['Retrieving NULL elements from the database (one time series ' ...
                'per database query). Please be patient...\n']);
case 'error'
    fprintf(1,['Retrieving error elements from the database (one time series ' ...
                'per database query). Please be patient...\n']);
end

% --------------------------------------------------------------------------
%% Retrieve the data from the database:
% --------------------------------------------------------------------------
iterationTimes = zeros(numTS,1); % Record the time taken for each iteration
didRetrieve = zeros(numTS,1);    % Keep track of whether data was retrieved at each iteration
for i = 1:numTS
    
    iterationTimer = tic; % Time each iteration using iterationTimer

    ts_id_now = ts_ids(i); % Range of ts_ids retrieved in this iteration
    
    % Start piecing together the mySQL SELECT command:
    switch retrieveWhatData
    case 'all' % retrieve all columns
        selectWhat = 'SELECT op_id, Output, QualityCode, CalculationTime FROM Results';
    case 'nocalctime' % Do not retrieve calculation time results
        selectWhat = 'SELECT op_id, Output, QualityCode FROM Results';
    case 'outputs' % just outputs
        selectWhat = 'SELECT op_id, Output FROM Results';
    case 'quality' % just quality codes
        selectWhat = 'SELECT op_id, QualityCode FROM Results';
    end
    
    baseString = sprintf('%s WHERE ts_id = %u AND op_id IN (%s)',selectWhat, ...
                                            ts_id_now,op_ids_string);
    
    % We could do a (kind of blind) retrieval, i.e., without retrieving op_ids safely
    % as long as for a given ts_id, the op_ids are in ascending order in the Results table.
    % This will be the case if time series and operations are added
    % using SQL_add because of the SORT BY specifier in SQL_add commands.
    % Otherwise op_ids should also be retrieved here, and used to sort
    % the other columns (i.e., outputs, quality codes, calculation times)
    
    switch retrieveWhatEntries
    case 'all'
        selectString = baseString;
    case 'null'
    	selectString = sprintf('%s AND QualityCode IS NULL',baseString);
    case 'error'
    	selectString = sprintf('%s AND QualityCode = 1',baseString);
    end
    
    % Do the retrieval
    % DatabaseTimer = tic;
	[qrc, emsg] = mysql_dbquery(dbc,selectString); % Retrieve data for this time series from the database
    % fprintf(1,'Database query for %u time series took %s\n',BundleSize,BF_thetime(toc(DatabaseTimer)));
    
    if ~isempty(emsg)
        error('Error retrieving outputs from %s.\n%s',dbname,emsg);
    end
    
    % Check results look ok:
    if isempty(qrc) % No data to retrieve
        % There are no entries in Results that match the requested conditions
        fprintf(1,'No data to retrieve for ts_id = %u\n',ts_id_now);
        % Leave local files (e.g., TS_DataMat, TS_Quality, TS_CalcTime as Inf)
        
    else
        % Entries need to be written to local matrices
        % Set didRetrieve = 1 for this iteration
        didRetrieve(i) = 1;
        
    	% Convert empty entries to NaNs
    	qrc(cellfun(@isempty,qrc)) = {NaN};
    
        % Put results from database into rows of local matrix
        if strcmp(retrieveWhatEntries,'all') % easy in this case
            % Assumes data is ordered by the op_id_string provided
            % This will be the case if all time series and operations
            % were added using SQL_add.
            % Otherwise we'll have nonsense happening...
            switch retrieveWhatData
            case 'all'
                TS_DataMat(i,:) = vertcat(qrc{:,2});
                TS_Quality(i,:) = vertcat(qrc{:,3});
                TS_CalcTime(i,:) = vertcat(qrc{:,4});
            case 'nocalctime'
                TS_DataMat(i,:) = vertcat(qrc{:,2});
                TS_Quality(i,:) = vertcat(qrc{:,3});
            case 'outputs'
                TS_DataMat(i,:) = vertcat(qrc{:,2});
            case 'quality'
                TS_Quality(i,:) = vertcat(qrc{:,2});
            end
        else
            % We retrieved a subset of the input op_ids
            % We have to match retrieved op_ids to local indicies
            iy = arrayfun(@(x)find(op_ids == x,1),vertcat(qrc{:,1}));
            % Now fill the corresponding entries in the local matrices:
            switch retrieveWhatData
            case 'all'
                TS_DataMat(i,iy) = vertcat(qrc{:,2});
                TS_Quality(i,iy) = vertcat(qrc{:,3});
                TS_CalcTime(i,iy) = vertcat(qrc{:,4});
            case 'nocalctime'
                TS_DataMat(i,iy) = vertcat(qrc{:,2});
                TS_Quality(i,iy) = vertcat(qrc{:,3});
            case 'outputs'
                TS_DataMat(i,iy) = vertcat(qrc{:,2});
            case 'quality'
                TS_Quality(i,iy) = vertcat(qrc{:,2});
            end
        end
    end
    
    % Note time taken for this iteration, and periodically display indication of time remaining
	iterationTimes(i) = toc(iterationTimer);
    if (i==1) % Give an initial indication of time after the first iteration
        fprintf(1,['Based on the first retrieval, this is taking ' ...
                'approximately %s per time series...\n'],BF_thetime(iterationTimes(1)));
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(iterationTimes(1)*(numTS-1)));
    elseif (mod(i,floor(numTS/10))==0) && i~=numTS % Tell us the time remaining 10 times across the total retrieval
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(mean(iterationTimes(1:i))*(numTS-i)));
	end
end

% --------------------------------------------------------------------------
%% Finished retrieving from the database!
% --------------------------------------------------------------------------
if any(didRetrieve)
    fprintf(1,'Retrieved data from %s over %u iterations in %s.\n',...
                            dbname,numTS,BF_thetime(sum(iterationTimes)));
else
    fprintf(1,['Over %u iterations, no data was retrieved from %s.\n' ...
                            'Not writing any data to file.\n'],numTS,dbname);
    SQL_closedatabase(dbc); return
end
	
if ismember(retrieveWhatEntries,{'null','error'})
    % We only want to keep rows and columns with (NaNs for 'null' or errors for 'error') in them...
    switch retrieveWhatEntries
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
	elseif sum(keepi) < numTS
		fprintf(1,'Cutting down from %u to %u time series\n',numTS,sum(keepi));
		ts_ids = ts_ids(keepi); numTS = length(ts_ids);
		ts_ids_string = BF_cat(ts_ids,',');
        switch retrieveWhatData
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
    elseif sum(keepi) < numOps
		fprintf(1,'Cutting down from %u to %u operations\n',numOps,sum(keepi));
		op_ids = op_ids(keepi); numOps = length(op_ids);
		op_ids_string = BF_cat(op_ids,',');
        switch retrieveWhatData
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


% ------------------------------------------------------------------------------
%% Fill Metadata
% ------------------------------------------------------------------------------

% 1. Retrieve Time Series Metadata
selectString = sprintf('SELECT FileName, Keywords, Length, Data FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string);
[tsinfo,emsg] = mysql_dbquery(dbc,selectString);
% Convert to a structure array, TimeSeries, containing metadata for all time series
tsinfo = [num2cell(ts_ids),tsinfo];
% Define inline functions to convert time-series data text to a vector of floats:
ScanCommas = @(x) textscan(x,'%f','Delimiter',',');
TakeFirstCell = @(x) x{1};
tsinfo(:,end) = cellfun(@(x) TakeFirstCell(ScanCommas(x)),tsinfo(:,end),'UniformOutput',0); % Do the conversion
TimeSeries = cell2struct(tsinfo',{'ID','FileName','Keywords','Length','Data'}); % Convert to structure array

% 2. Retrieve Operation Metadata
selectString = sprintf('SELECT OpName, Keywords, Code, mop_id FROM Operations WHERE op_id IN (%s)',op_ids_string);
[opinfo,emsg] = mysql_dbquery(dbc,selectString);
opinfo = [num2cell(op_ids), opinfo]; % add op_ids
Operations = cell2struct(opinfo',{'ID','Name','Keywords','CodeString','MasterID'});

% Check that no operations have bad links to master operations:
if any(isnan([Operations.MasterID]))
    fisBad = find(isnan([Operations.MasterID]));
    for i = 1:length(fisBad)
        fprintf(1,'Bad link (no master match): %s -- %s\n',Operations(fisBad(i)).Name,Operations(fisBad(i)).CodeString);
    end
    error('Bad links of %u operations to non-existent master operations',length(fisBad));
end

% 3. Retrieve Master Operation Metadata
% (i) Which masters are implicated?
selectString = ['SELECT mop_id, MasterLabel, MasterCode FROM MasterOperations WHERE mop_id IN ' ...
        				'(' BF_cat(unique([Operations.MasterID]),',') ')'];
[masterinfo,emsg] = mysql_dbquery(dbc,selectString);
if ~isempty(emsg)
    fprintf(1,'Error retrieving Master information using:\n%s\n',selectString);
    disp(emsg);
    error('Master information could not be retrieved')
else
    MasterOperations = cell2struct(masterinfo',{'ID','Label','Code'});
end

% Close database connection
SQL_closedatabase(dbc)

% ------------------------------------------------------------------------------
%% Save to HCTSA_loc.mat
% ------------------------------------------------------------------------------
fprintf(1,'Saving local versions of the data to HCTSA_loc.mat...');
save('HCTSA_loc.mat','TimeSeries','Operations','MasterOperations','-v7.3');
switch retrieveWhatData
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
didWrite = 1; % Tag to say that write to HCTSA_loc.mat file is successful

% ------------------------------------------------------------------------------
%% If retrieved quality labels, display how many entries need to be calculated
% ------------------------------------------------------------------------------
if strcmp(retrieveWhatData,'outputs')
    fprintf(1,'You have the outputs, but you don''t know which are good or not without the quality labels...\n');
else
    toCalculate = sum(isnan(TS_Quality(:)) | TS_Quality(:)==1);
    fprintf(1,'There are %u entries (=%4.2f%%) to calculate in the data matrix (%ux%u).\n', ...
                                        toCalculate,toCalculate/numOps/numTS*100,numTS,numOps);
end

end