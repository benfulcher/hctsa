function didWrite = SQL_retrieve(ts_ids,op_ids,retrieveWhatEntries,retrieveWhatData)
% SQL_retrieve 		Retrieve data from the mySQL database
%
% This function retreives data from the mySQL database for subsequent analysis
% in Matlab. It takes as input a set of constraints on the time series and
% operations to include, and outputs the relevant subsection of the data matrix
% and associated metadata to HCTSA.mat
%
%---INPUTS:
% ts_ids: a vector of ts_ids to retrieve from the mySQL database.
% op_ids: a vector of op_ids to retrieve from the mySQL database.
% retrieveWhatEntries: can be one of the following options:
%        (i) 'null': Retrieve null entries from the database (e.g., for
%                    calculating)
%        (ii) 'all': Retrieve all entries from the database (e.g., for
%                    analyzing)
%        (iii) 'error': Retrieve previous errors stored in the database (e.g.,
%                       for re-evaluating previous errors)
% retrieveWhatData: can be one of the following options:
%       (i) 'all': Retrieves Data, CalcTimes, and Quality labels
%       (ii) 'nocalctime': Retrieves calculated values and quality labels
%       (iii) 'outputs': Just the calculated values
%       (iv) 'quality': Just the quality labels
%
%---OUTPUT:
% didWrite [opt] is 1 if HCTSA.mat was written
%
% Other outputs are saved to the file HCTSA.mat:
%--TS_DataMat, contains the data
%--TS_Quality, contains the quality codes
%--TS_CalcTime, contains the calculation times [if retrieveWhatData = 'all']
%--TimeSeries, contains metadata about the time series, and the time-series data
%--Operations, contains metadata about the operations
%--MasterOperations, contains metadata about the implicatedmaster operations

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Until it actually writes something, set the function output, didWrite = 0
didWrite = 0;

% --------------------------------------------------------------------------
%% Check inputs and set defaults
% --------------------------------------------------------------------------
if nargin < 2
	error('You must provide at least two inputs!');
end
if (ischar(ts_ids) && ~strcmp(ts_ids,'all')) || (ischar(op_ids) && ~strcmp(op_ids,'all'))
	error('Must specify ''all'' when retrieving IDs, or specify a custom set as a vector');
end

% retrieveWhatEntries
if nargin < 3 || isempty(retrieveWhatEntries)
	retrieveWhatEntries = 'all'; % default: retrieve full sets of things, not just the empty entries in the database
    fprintf(1,'Retrieving ALL elements from the database in the specified ts_id and op_id ranges by default.\n');
end
retrieveWhatEntriesCanBe = {'null','all','error'};
if ~ischar(retrieveWhatEntries) || ~ismember(retrieveWhatEntries,retrieveWhatEntriesCanBe)
    error(['The third input to SQL_retrieve must specify what type of entries to retrieve from ' ...
                'the database, one of the following: %s'],BF_cat(retrieveWhatEntriesCanBe))
end

% retrieveWhatData
if nargin < 4 || isempty(retrieveWhatData)
    retrieveWhatData = 'nocalctime'; % default: just retrieve data and quality labels
end
retrieveWhatDataCanBe = {'all','nocalctime','outputs','quality'};
if ~ischar(retrieveWhatData) || ~ismember(retrieveWhatData,retrieveWhatDataCanBe)
    error(['The fourth input to SQL_retrieve must specify what data to retrieve from the database, ' ...
                'one of the following: %s'],BF_cat(retrieveWhatEntriesCanBe))
end

% ------------------------------------------------------------------------------
% First check whether you're about to overwrite an existing file
% (actually you probably don't want to do this when you're looping in a script)
% ------------------------------------------------------------------------------
% if exist('./HCTSA.mat','file')
%     reply = input(['Warning: HCTSA.mat already exists -- if you continue, this ' ...
%                     'file will be overwritten.\n[press any key to continue]'],'s');
% end

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------

% Process ts_ids if provided a list (i.e., not 'all')
if ~(ischar(ts_ids) && strcmp(ts_ids,'all'))
	[ts_ids,ts_ids_string,numTS] = idProcessing(ts_ids);
	if numTS == 0
		error('Oops. There''s nothing to do! No time series to retrieve!');
	end
end

if ~(ischar(op_ids) && strcmp(op_ids,'all'))
	[op_ids,op_ids_string,numOps] = idProcessing(op_ids);
	if numOps == 0
		error('Oops. There''s nothing to do! No operations to retrieve!');
	end
end

% ------------------------------------------------------------------------------
% Open a database connection
% ------------------------------------------------------------------------------
% (attempt to use the Matlab database toolbox, which is faster for retrievals)
[dbc, dbname] = SQL_opendatabase('',false,true);

% ------------------------------------------------------------------------------
% Refine the set of time series and operations to those that actually exist in
% the database
% ------------------------------------------------------------------------------

% TimeSeries IDs that exist in the database:
if ischar(ts_ids) && strcmp(ts_ids,'all')
	tsids_db = mysql_dbquery(dbc,'SELECT ts_id FROM TimeSeries');
	tsids_db = vertcat(tsids_db{:}); % Will always be sorted in ascending order
else
	tsids_db = mysql_dbquery(dbc,sprintf('SELECT ts_id FROM TimeSeries WHERE ts_id IN (%s)',ts_ids_string));
	tsids_db = vertcat(tsids_db{:});
	% There are fewer time series in the database than requested
	if length(tsids_db) < numTS
	    if isempty(tsids_db) % Now there are no time series to retrieve
	        fprintf(1,'None of the %u specified time series exist in ''%s''\n',numTS,dbname);
	        SQL_closedatabase(dbc); return % Close the database connection before returning
	    end
	    fprintf(1,['%u specified time series do not exist in ''%s'', retrieving' ...
	                    ' the remaining %u\n'],numTS-length(tsids_db),dbname,length(tsids_db));
	end
end
ts_ids_string = BF_cat(tsids_db,',');
numTS = length(tsids_db);

% Operation IDs that exist in the database:
if ischar(op_ids) && strcmp(op_ids,'all')
	opids_db = mysql_dbquery(dbc,'SELECT op_id FROM Operations');
	opids_db = vertcat(opids_db{:});
else
	opids_db = mysql_dbquery(dbc,sprintf('SELECT op_id FROM Operations WHERE op_id IN (%s)',op_ids_string));
	opids_db = vertcat(opids_db{:});
	% There are fewer operations in the database than requested:
	if length(opids_db) < numOps
	    if isempty(opids_db) % Now there are no operations to retrieve
	        fprintf(1,'None of the %u specified operations exist in ''%s''\n',numOps,dbname);
	        SQL_closedatabase(dbc); return % Close the database connection before returning
	    end
	    fprintf(1,['%u specified operations do not exist in ''%s'', retrieving' ...
	                    ' the remaining %u\n'],numOps-length(opids_db),dbname,length(opids_db));
	end
end
op_ids_string = BF_cat(opids_db,',');
numOps = length(opids_db);

%-------------------------------------------------------------------------------
% Tell me about it:
fprintf(1,'We have %u time series and %u operations to retrieve from %s.\n',...
						numTS,numOps,dbname);
fprintf(1,['Filling and saving to local Matlab file HCTSA.mat from ' ...
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

didRetrieve = zeros(numTS,1); % Keep track of whether data was retrieved at each iteration
retrievalTimer = tic; % Time the process using retrievalTimer

for i = 1:numTS

    ts_id_now = tsids_db(i); % Range of ts_ids retrieved in this iteration

	if ischar(op_ids) && strcmp(op_ids,'all')
		% All op_ids
		baseString = sprintf('%s WHERE ts_id = %u',selectWhat,ts_id_now);
	else
		% Custom set of op_ids
	    baseString = sprintf('%s WHERE ts_id = %u AND op_id IN (%s)',selectWhat, ...
                                            ts_id_now,op_ids_string);
	end

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
	[qrc, emsg] = mysql_dbquery(dbc,selectString); % Retrieve data for this time series from the database

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
            iy = arrayfun(@(x)find(opids_db == x,1),vertcat(qrc{:,1}));
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

    % Periodically display indication of time remaining
    if (i==50) && (i < numTS/10) % Give an initial indication of time after the first 50 iterations
        fprintf(1,['Based on the first 50 retrievals, this is taking ' ...
                'approximately %s per time series...\n'],BF_thetime(toc(retrievalTimer)/i));
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(toc(retrievalTimer)/i*(numTS-i)));
    elseif (mod(i,floor(numTS/10))==0) && i~=numTS % Tell us the time remaining 10 times across the total retrieval
		fprintf(1,'Approximately %s remaining...\n',BF_thetime(toc(retrievalTimer)/i*(numTS-i)));
	end
end

% --------------------------------------------------------------------------
%% Finished retrieving from the database!
% --------------------------------------------------------------------------
if any(didRetrieve)
    fprintf(1,'Retrieved data from %s over %u iterations in %s.\n',...
                            dbname,numTS,BF_thetime(toc(retrievalTimer)));
else
    fprintf(1,['Over %u iterations, no data was retrieved from %s (%s).\n' ...
                            'Not writing any data to file.\n'],...
							numTS,dbname,BF_thetime(toc(retrievalTimer)));
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
		tsids_db = tsids_db(keepi);
        numTS = length(tsids_db);
		ts_ids_string = BF_cat(tsids_db,',');
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
		opids_db = opids_db(keepi);
		numOps = length(opids_db);
		op_ids_string = BF_cat(opids_db,',');
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
if ischar(ts_ids) && strcmp(ts_ids,'all')
	selectString = 'SELECT Name, Keywords, Length, Data FROM TimeSeries';
else
	selectString = sprintf(['SELECT Name, Keywords, Length, Data FROM ',...
								'TimeSeries WHERE ts_id IN (%s)'],ts_ids_string);
end
[tsinfo,emsg] = mysql_dbquery(dbc,selectString);
if ~isempty(emsg)
    error('Error retrieving time-series metadata from from %s',dbname);
end
% Convert to a structure array, TimeSeries, containing metadata for all time series
tsinfo = [num2cell(tsids_db),tsinfo];
% Define inline functions to convert time-series data text to a vector of floats:
scanCommas = @(x) textscan(x,'%f','Delimiter',',');
takeFirstCell = @(x) x{1};
tsinfo(:,end) = cellfun(@(x) takeFirstCell(scanCommas(x)),tsinfo(:,end),'UniformOutput',0); % Do the conversion
TimeSeries = cell2struct(tsinfo',{'ID','Name','Keywords','Length','Data'}); % Convert to structure array


% 2. Retrieve Operation Metadata
% (even if specify 'all', can be fewer for 'null','error' cases, where you restrict
% the operations that are actually retrieved to the local file)
% [would still probably be faster to retrieve all above, and then subset the info using keepi]
selectString = sprintf('SELECT Name, Keywords, Code, mop_id FROM Operations WHERE op_id IN (%s)',op_ids_string);
opinfo = mysql_dbquery(dbc,selectString);
opinfo = [num2cell(opids_db), opinfo]; % add op_ids
Operations = cell2struct(opinfo',{'ID','Name','Keywords','CodeString','MasterID'});


% Check that no operations have bad links to master operations:
if any(isnan([Operations.MasterID]))
    fisBad = find(isnan([Operations.MasterID]));
    for i = 1:length(fisBad)
        fprintf(1,'Bad link (no master match): %s -- %s\n',...
					Operations(fisBad(i)).Name,Operations(fisBad(i)).CodeString);
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
%% Save to HCTSA.mat
% ------------------------------------------------------------------------------
outputFileName = 'HCTSA.mat';
fprintf(1,'Saving local versions of the data to %s...',outputFileName);
saveTimer = tic;
fromDatabase = 1; % mark that we retrieved this data from the mySQL database
save('HCTSA.mat','TimeSeries','Operations','MasterOperations','fromDatabase','-v7.3');
switch retrieveWhatData
case 'all'
    % Add outputs, quality labels, and calculation times
    save(outputFileName,'TS_DataMat','TS_Quality','TS_CalcTime','-append')
case 'nocalctime'
    % Add outputs and quality labels
    save(outputFileName,'TS_DataMat','TS_Quality','-append')
case 'outputs'
    % Add outputs
    save(outputFileName,'TS_DataMat','-append')
case 'quality'
    % Add quality labels
    save(outputFileName,'TS_Quality','-append')
end

fprintf(1,' Done in %s.\n',BF_thetime(toc(saveTimer)));
clear saveTimer % stop timing

didWrite = 1; % Tag to say that write to HCTSA.mat file is successful

% ------------------------------------------------------------------------------
%% If retrieved quality labels, display how many entries need to be calculated
% ------------------------------------------------------------------------------
if strcmp(retrieveWhatData,'outputs')
    fprintf(1,'You have the outputs, but you don''t know which are good or not without the quality labels...\n');
else
    toCalculate = sum(isnan(TS_Quality(:)) | TS_Quality(:)==1);
    fprintf(1,'There are %u entries (%4.2f%%) to calculate in the %u x %u data matrix.\n', ...
                                        toCalculate,toCalculate/numOps/numTS*100,numTS,numOps);
end


%-------------------------------------------------------------------------------
function [idList,idString,numIDs] = idProcessing(idList)
	% Make sure id lists are column vectors:
	if size(idList,2) > size(idList,1), idList = idList'; end
	% Sort ids ascending:
	idList = sort(idList,'ascend');
	% Write a comma-delimited string of ids:
	idString = BF_cat(idList,',');
	% Count the number:
	numIDs = length(idList);
end

end
