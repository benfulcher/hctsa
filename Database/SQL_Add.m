function theTable = SQL_Add(addWhat,inputFile,forDatabase,beVocal)
% SQL_Add   Interpret a structured input file of time series, operations,
%           or master operations.
%
% By default adds the results to a linked mySQL database.
%
%---INPUTS:
% addWhat: 'mops' (for master operations), 'ops' (for operations), or 'ts'
%             (for time series)
% inputFile: the filename of the tab-delimited textfile to be read in [default
%            = INP_ts.txt or INP_ops.txt or INP_mops.txt]
%            The input file should be formatted with whitespace as a delimiter
%            between the entries to import.
% forDatabase: if true (default) write the results of interpreting the input file
%               to the default linked mySQL database for hctsa.
% beVocal: if true (default) gives user feedback on the input process.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------

% addWhat: 'ts', 'mops', or 'ops'
if nargin < 1 || isempty(addWhat) || ~ismember(addWhat,{'ops','ts','mops'})
    error(['Error setting first input argument -- should be ''ts'' for TimeSeries ' ...
                ', ''ops'' for Operations, or ''mops'' for Master Operations']);
end

% inputFile
if nargin < 2 || isempty(inputFile)
    % Default filenames:
    switch inputFile
    case 'ts'
        inputFile = 'INP_ts.txt';
    case 'ops'
        inputFile = 'INP_ops.txt';
    case 'mops'
        inputFile = 'INP_mops.txt';
    end
end
if ~exist(inputFile,'file')
    error('Unknown file ''%s''',inputFile);
end

if nargin < 3 || isempty(forDatabase)
    forDatabase = false;
end
if nargin < 4
    beVocal = true; % Give user feedback by default:
end

% ------------------------------------------------------------------------------
% Display welcome message:
% ------------------------------------------------------------------------------
if beVocal
    fprintf(1,'Using input file: %s\n',inputFile);
end
ticker = tic;

% ------------------------------------------------------------------------------
%% Open database connection
% ------------------------------------------------------------------------------
if forDatabase
    [dbc,databaseName] = SQL_OpenDatabase();
end

% ------------------------------------------------------------------------------
% Define strings to unify the different strands of code for time series /
% operations
% ------------------------------------------------------------------------------
switch addWhat
    case 'ts'
        % Check that the time series table exists in the database:
        theWhat = 'time series';
        theid = 'ts_id';
        thekid = 'tskw_id';
        theTable = 'TimeSeries';
        theKeywordTable = 'TimeSeriesKeywords';
        theRelTable = 'TsKeywordsRelate';
        maxL = 50000; % the longest time series length accepted in the database
    case 'ops'
        theWhat = 'operations';
        theid = 'op_id';
        thekid = 'opkw_id';
        theTable = 'Operations';
        theKeywordTable = 'OperationKeywords';
        theRelTable = 'OpKeywordsRelate';
    case 'mops'
        theWhat = 'master operations';
        theid = 'mop_id';
        theTable = 'MasterOperations';
end

% ------------------------------------------------------------------------------
% Check that the table exists in the datbase
% ------------------------------------------------------------------------------
if forDatabase
    existString = ['SHOW TABLES LIKE ''' theTable ''''];
    [output,emsg] = mysql_dbquery(dbc,existString);
    if isempty(output) % Table doesn't exist
        error(['Table %s doesn''t exist in the database %s.\n' ...
                    'Use: (1) install.m to set the database system up from scratch,\n' ...
                    '(2) SQL_CreateAllTables to create empty tables that can later be filled with custom libraries of operations,\n' ...
                    '(3) SQL_Reset to drop all tables in the database, and repopulate all tables with the default library of operations.'],...
                    theTable,databaseName);
    end
end

% ------------------------------------------------------------------------------
%% Open and read the input file
% ------------------------------------------------------------------------------

% Determine if input is a .mat file or an input text file based on filename extension:
if strcmp(inputFile(end-3:end),'.mat');
    isMatFile = true;
else
    isMatFile = false;
end

if ~isMatFile
    % Specified a input text file
    fid = fopen(inputFile);
    if (fid==-1)
        error('Could not load the specified input text file ''%s''',inputFile)
    end
    switch addWhat
    case 'ts' % Read the time series input file:
        if beVocal
            fprintf(1,['Need to format %s (Time Series input file) as: Name ' ...
                                                                'Keywords\n'],inputFile);
            fprintf(1,'Assuming no header line\n');
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n');
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n');
        end
        formatSpec = '%s%s';
    case 'ops' % Read the operations input file:
        if beVocal
            fprintf(1,['Need to format %s (Operations input file) as: OperationCode ' ...
                                            'OperationName OperationKeywords\n'],inputFile);
            fprintf(1,'Assuming no header lines\n');
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n');
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n');
        end
        formatSpec = '%s%s%s';
    case 'mops' % Read the master operations input file:
        if beVocal
            fprintf(1,'Need to format %s (Master Operations input file) as: MasterCode MasterLabel\n',inputFile);
            fprintf(1,'Assuming no header lines\n');
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n');
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n');
        end
        formatSpec = '%s%s';
    end
    dataIn = textscan(fid,formatSpec,'CommentStyle','#','EndOfLine','\r\n','CollectOutput',true);
    fclose(fid);

    % ------------------------------------------------------------------------------
    % Show the user what's been imported:
    % ------------------------------------------------------------------------------
    dataIn = dataIn{1}; % Collect one big matrix of cells
    numItems = size(dataIn,1); % Number of items in the input file

    if numItems == 0
        error('The input file ''%s'' seems to be empty??',inputFile)
    end

    if beVocal
        fprintf(1,'Found %u %s in %s, I think. Take a look:\n\n',numItems,theWhat,inputFile);
        switch addWhat
        case 'ts'
            fprintf(1,'%s\t%s\n','-Name-','-Keywords-');
            fprint_ts = @(x)fprintf('%s\t%s\n',dataIn{x,1},dataIn{x,2});
        case 'ops'
            fprintf(1,'%s\t%s\t%s\n','-Operation Name-','-Master Label-','-Operation Keywords-');
            fprint_ops = @(x) fprintf('%s\t%s\t%s\n',dataIn{x,1},dataIn{x,2},dataIn{x,3});
        case 'mops'
            fprintf(1,'%s\t%s\n','-Master Code-','-Master Label-');
            fprint_mops = @(x) fprintf('%s\t%s\n',dataIn{x,1},dataIn{x,2});
        end

        for i = 1:min(3,numItems)
            switch addWhat
            case 'ts', fprint_ts(i);
            case 'ops', fprint_ops(i);
            case 'mops', fprint_mops(i);
            end
        end

        if numItems > 3
            fprintf(1,'..................(%u).....................\n',max(numItems-6,0));
            for i = max(numItems-2,4):numItems
                switch addWhat
                case 'ts', fprint_ts(i);
                case 'ops', fprint_ops(i);
                case 'mops', fprint_mops(i);
                end
            end
        end

        fprintf(1,['\nHow does it look? Make sure the metadata matches the headings\n']);

        % Ask the question:
        if strcmp(addWhat,'ts')
            if forDatabase
                reply = input(['If we go on, we will attempt to read all time series ' ...
                            'from file and add all ' ...
                            'data to the database.\n<<<Type ''y'' to continue...>>> '],'s');
            else
                reply = input(['If we go on, we will attempt to read all time series ' ...
                            'from file and add all ' ...
                            'data to HCTSA.mat\n<<<Type ''y'' to continue...>>> '],'s');
            end
        else
            if forDatabase
                reply = input(['If we go on, we will attempt to add all ' ...
                            'data to the database.\n<<<Type ''y'' to continue...>>> '],'s');
            else
                reply = input(['If we go on, we will attempt to add all ' ...
                            'data to HCTSA.mat\n<<<Type ''y'' to continue...>>> '],'s');
            end
        end

        if ~strcmp(reply,'y')
            fprintf(1,'I didn''t think so. Come back later...\n');
            return
        end
    end

    if beVocal
        fprintf(1,'%s read.\n',inputFile);
    end

else
    % .mat file input is only allowed when importing time series
    if ~strcmp(addWhat,'ts')
        error(['.mat file input type only supported for importing time series. ' ...
                'Please specify a text file for importing master operations or operations.']);
    end

    % Load the 3 cells specifying the data:
    if ~exist(inputFile,'file')
        error('Could not load specified input .mat file: %s',inputFile);
    end
    inputData = load(inputFile,'timeSeriesData','labels','keywords');

    % Check that they're all as they should be:
    if ~isfield(inputData,'timeSeriesData') || ~isfield(inputData,'labels') || ~isfield(inputData,'keywords') ...
                || ~(iscell(inputData.timeSeriesData) || isnumeric(inputData.timeSeriesData)) ...
                || ~iscell(inputData.labels) || ~iscell(inputData.keywords)
        error(['Incorrectly formatted input file, %s\nExpecting input file to ' ...
                    'contain: ''timeSeriesData'' (cell or matrix), ' ...
                    '''labels'' (cell), and ''keywords'' (cell).'],inputFile);
    end

    % Get number of time series (as numItems):
    if iscell(inputData.timeSeriesData)
        numItems = length(inputData.timeSeriesData);
    else
        % Specified a matrix of time-series data
        numItems = size(inputData.timeSeriesData,1);
    end
    if beVocal
        fprintf(1,['We have %u time series (rows) in your data matrix, timeSeriesData ' ...
                                    '(loaded from %s).\n'],numItems,inputFile);
        if forDatabase
            fprintf(1,['Will store time-series data from matlab file in the database to an' ...
                        ' accuracy of 6 significant figures...\n']);
        end
    end

    % Check sizes match:
    if length(inputData.labels) ~= numItems || length(inputData.keywords) ~= numItems
        fprintf(1,'%u time series, %u labels, %u keywords in %s.\n',numItems,...
                    length(inputData.labels),length(inputData.keywords),inputFile);
        error('All cells in the input file %s must be the same length.',inputFile);
    end

    % Plot some to screen
    if beVocal
        plotNum = min(5,numItems); % display as many as the first 5 time series to demonstrate
        figure('color','w','WindowStyle','docked');
        for j = 1:plotNum
            subplot(plotNum,1,j);
            if iscell(inputData.timeSeriesData)
                plot(inputData.timeSeriesData{j},'k');
            else
                plot(inputData.timeSeriesData(j,:),'k');
            end
            title(sprintf('[%u/%u] %s (%s)',j,numItems,inputData.labels{j},inputData.keywords{j}),'interpreter','none')
            if j==plotNum
                xlabel('Time (samples)')
            end
        end

        if forDatabase
            reply = input(sprintf(['Does this look ok for the first %u time series?\nIf we continue, ' ...
                'we will attempt to add all %u time series in the input file' ...
                ' to the database.\n<<<Type ''y'' to continue...>>> '],plotNum,numItems),'s');
        else
            reply = input(sprintf(['Does this look ok for the first %u time series?\n' ...
                                    '<<<Type ''y'' to continue...>>> '],plotNum),'s');
        end
        if ~strcmp(reply,'y')
            fprintf(1,'I didn''t think so. Come back later...\n');
            return
        end
        close
    end

    % Ok, so we have inputData.timeSeriesData, inputData.labels, and inputData.keywords
end

% ------------------------------------------------------------------------------
% Construct a table of metadta for the time series / operations /
% master operations and fill a cell, toAdd, containing mySQL INSERT commands for
% each item in the input file.
% ------------------------------------------------------------------------------
if forDatabase
    % Define an inline function to add escape strings to format mySQL queries
    if forDatabase
        esc = @RA_sqlescapestring;
        % Initialize:
        toAdd = cell(numItems,1);
    end
    if beVocal
        fprintf(1,'Preparing mySQL statements to add %u %s to the database %s...',...
                    numItems,theWhat,databaseName);
    end
end

% Prepare local IDs:
ID = (1:numItems)';

switch addWhat
case 'ts' % Prepare toAdd cell for time series
    if beVocal
        figure('color','w','WindowStyle','docked');
    end
    % Record whether data was added or not (not too long and no missing values)
    wasGood = false(1,numItems);

    % Initialize:
    Name = cell(numItems,1);
    Keywords = cell(numItems,1);
    Length = zeros(numItems,1);
    Data = cell(numItems,1);
    for j = 1:numItems
        % Assign filename and keywords strings to this time series, and load it
        if isMatFile
            Name{j} = inputData.labels{j};
            Keywords{j} = inputData.keywords{j}; % Take out inverted commas from keywords lists
            if iscell(inputData.timeSeriesData)
                Data{j} = inputData.timeSeriesData{j};
                if size(Data{j},2) > size(Data{j},1)
                    fprintf(1,'Transposing time series\n');
                    Data{j} = Data{j}';
                end
                if size(Data{j},2) ~= 1
                    error('Multivariate time-series input? Each element of timeSeriesData must be univariate');
                end
            else
                Data{j} = inputData.timeSeriesData(j,:)'; % column vector
            end
        else
            Name{j} = dataIn{j,1};
            Keywords{j} = regexprep(dataIn{j,2},'\"',''); % Take out inverted commas from keywords lists
            % Read the time series from its filename:
            try
                Data{j} = dlmread(Name{j});
            catch emsg
                fprintf(1,'%s\n',emsg.message);
                error(['\nCould not read the data file for ''%s''.',...
                        ' Check that it''s in the path.'],Name{j})
            end
        end

        % Check that label is not empty:
        if isempty(Name{j})
            warning(['\n[%u/%u] This time series is assigned an empty label' ...
                        ' and will not be added...'],j,numItems)
            beep; continue
        end

        % Assign the length of the time series
        Length(j) = length(Data{j});

        % TEST 1: Is the time series longer than the maximum allowed value?
        if Length(j) > maxL
            warning(['\n[%u/%u]%s contains %u samples, this framework can efficiently ' ...
                'deal with time series up to %u samples\n',...
                '[Note that this maximum length can be modified in SQL_Add]\n',...
                'Skipping this time series...'],...
                j,numItems,Name{j},Length(j),maxL)
            beep; continue
        end
        % TEST 2: Does the time series contains any NaN of Inf values?
        if any(isnan(Data{j})) || any(~isfinite(Data{j}))
            warning(['\n[%u/%u] The time series: ''%s'' contains special values' ...
                        ' (e.g., NaN or Inf)...\nThis time series will not be added...'], ...
                        j,numItems,Name{j})
            beep; continue
        end

        % Passed both tests! Assign wasGood = true
        wasGood(j) = true;

        % If storing in a database, need to reassign the time-series data as text
        % (which will be stored as singles in the case of a .mat file data)
        % If making for a local mat file, keep the same precision as input.
        if forDatabase
            % Want time-series data as a comma-delimited string:
            if isMatFile
                xText = sprintf('%.6g,',Data{j}); % keeps 6 figures of accuracy
                xText = xText(1:end-1); % remove trailing comma
            else
                % Read in the time series from file as strings using textscan
                fid = fopen(Name{j});
                timeSeriesData_text = textscan(fid,'%s');
                fclose(fid);
                % Turn the time series into a comma-delimited string
                timeSeriesData_text = timeSeriesData_text{1}; % cell of time series values as strings
                xText = ''; % make a comma-delimited text version of the time series as xText
                for k = 1:length(timeSeriesData_text)
                    xText = [xText,',',timeSeriesData_text{k}];
                end
                xText = xText(2:end);
            end
            Data{j} = xText;
            % Prepare the data to be added to the database in an INSERT command:
            toAdd{j} = sprintf('(''%s'',''%s'',%u,''%s'')',...
                    esc(Name{j}),esc(Keywords{j}),Length{j},Data{j});
        end

        if beVocal
            % Plot the time series
            numSubplots = min(numItems,4);
            ax = subplot(numSubplots,1,mod(j-1,numSubplots)+1);
            plot(Data{j},'-k');
            ax.XLim = [1,Length(j)];
            titleText = sprintf('[%u/%u] %s (%u), keywords = %s',j,numItems,...
                            Name{j},Length(j),Keywords{j});
            title(titleText,'interpreter','none');
            fprintf(1,'\n%s --- loaded successfully.',titleText);
            pause(5e-3); % Wait 5ms to show the plotted time series!
        end
    end

    %-------------------------------------------------------------------------------
    % All data in -> TimeSeries table
    TimeSeries = table(ID,Name,Keywords,Length,Data);

    % Check for duplicates on the name field:
    if length(unique(TimeSeries.Name)) < height(TimeSeries)
        error('Input file contains duplicate names.');
    end

    % Check whether at least some passed the quality checks
    if ~any(wasGood)
        fprintf(1,'NONE of the %u time series in the input file passed quality checks.\n',...
                                length(wasGood));
        return
    end

    if beVocal
        % Tell me about it:
        if forDatabase
            textShow = ', ready to be uploaded to the mySQL database';
        else
            textShow = '';
        end
        fprintf(1,'\nAll time-series data loaded; %u/%u passed quality tests%s.\n',...
                        sum(wasGood),length(wasGood),textShow);
        if any(~wasGood)
            input(sprintf('[List %u time series that failed... (press any key)]',sum(~wasGood)),'s');
            iNoGood = find(~wasGood);
            for i = 1:length(iNoGood)
                if ~isempty(TimeSeries.Length(iNoGood(i)))
                    lengthText = sprintf(', N = %u',TimeSeries.Length(iNoGood(i)));
                else
                    lengthText = '';
                end
                fprintf(1,'*NOT UPLOADING:* [%u] %s (%s)%s\n',iNoGood(i),...
                    TimeSeries.Name{iNoGood(i)},TimeSeries.Keywords{iNoGood(i)},...
                    lengthText);
            end
            input(sprintf('[press any key to continue to add the remaining %u time series]',...
                                                        sum(wasGood)),'s');
        end
    end

case 'mops'
    % Prepare metadata for master operations
    Code = cell(numItems,1);
    Label = cell(numItems,1);
    for j = 1:numItems
        Code{j} = dataIn{j,1};
        Label{j} = dataIn{j,2};
        if forDatabase
            toAdd{j} = sprintf('(''%s'', ''%s'')',esc(Label{j}),esc(Code{j}));
        end
    end
    % Convert to table of metadata:
    MasterOperations = table(ID,Label,Code);
    if beVocal, fprintf(1,' Done.\n'); end

case 'ops'
    % Prepare metadata table for operations
    CodeString = cell(numItems,1);
    Name = cell(numItems,1);
    Keywords = cell(numItems,1);
    Label = cell(numItems,1);
    for j = 1:numItems
        CodeString{j} = dataIn{j,1};
        Name{j} = dataIn{j,2};
        Keywords{j} = dataIn{j,3};
        Label{j} = strtok(CodeString{j},'.');
        if forDatabase
            toAdd{j} = sprintf('(''%s'', ''%s'',''%s'',''%s'')',...
                    esc(Name{j}),esc(CodeString{j}),...
                    esc(Label{j}),esc(Keywords{j}));
        end
    end
    Operations = table(ID,Name,Label,Keywords,CodeString);
    if beVocal, fprintf(1,' Done.\n'); end

    % Check for duplicates in the input file:
    [uniqueOpNames,ia] = unique(Operations.Name);
    if length(uniqueOpNames) < height(Operations)
        warning(['Input file contains %u duplicate entries, which are being removed.\n' ...
                        'Inputting %u -> %u operations...'], ...
            height(Operations)-length(uniqueOpNames),...
            height(Operations),length(uniqueOpNames));
        % Check for duplicates:
        nameTable = tabulate(Operations.Name);
        isDuplicate = find([nameTable{:,2}] > 1);
        fprintf(1,'E.g., %s is a duplicate\n',nameTable{isDuplicate(1),1});
        % Only keep the unique ones:
        Operations = Operations(ia,:);
        numItems = height(Operations);
        fprintf(1,'We now have %u operations to input...\n',numItems);
        if forDatabase
            toAdd = toAdd(ia);
        end
    end
end

% ------------------------------------------------------------------------------
% If not writing to a database, then we're done
% ------------------------------------------------------------------------------
switch addWhat
case 'ts'
    theTable = TimeSeries(wasGood,:);
    % Note that IDs are according to the original input file
case 'mops'
    theTable = MasterOperations;
case 'ops'
    theTable = Operations;
end
if ~forDatabase
    fprintf(1,'Returning a table with metadata for %u %s.\n',height(theTable),theWhat);
    return
end

%-------------------------------------------------------------------------------
% From here on is for mySQL data import:
% ------------------------------------------------------------------------------
% Check for duplicates:
if beVocal
    fprintf(1,'Checking for duplicates already in the database... ');
end
switch addWhat
case 'ts'
    existingNames = mysql_dbquery(dbc,sprintf('SELECT Name FROM TimeSeries'));
    isDuplicate = ismember(theTable.Name,existingNames);
case 'ops'
    existingOperationNames = mysql_dbquery(dbc,'SELECT Name FROM Operations');
    isDuplicate = ismember(theTable.Name,existingOperationNames);
case 'mops'
    existingMasterLabels = mysql_dbquery(dbc,'SELECT MasterLabel FROM MasterOperations');
    isDuplicate = ismember(theTable.Label,existingMasterLabels);
end
if beVocal, fprintf(1,'Done.\n'); end

% ------------------------------------------------------------------------------
% Tell the user about duplicates
% ------------------------------------------------------------------------------
if all(isDuplicate)
    fprintf(1,'All %u %s from %s already exist in %s---no new %s to add!\n',...
                    numItems,theWhat,inputFile,databaseName,theWhat);
    return
elseif sum(isDuplicate) > 0
    if beVocal
        fprintf(1,'I found %u duplicate %s already in the database %s!\n',...
                        sum(isDuplicate),theWhat,databaseName);
        fprintf(1,'There are %u new %s to add to %s.\n',...
                        sum(~isDuplicate),theWhat,databaseName);
    end
end

% ------------------------------------------------------------------------------
% Incorporate isDuplicate with wasGood
% ------------------------------------------------------------------------------
if strcmp(addWhat,'ts')
    isBad = (isDuplicate | ~wasGood);
    if all(isBad)
        fprintf(1,['All input time series either already exist in the database, ' ...
                    'or did not pass the required quality checks. Exiting.']);
        return
    end
else
    isBad = isDuplicate;
end

% ------------------------------------------------------------------------------
%% Select the maximum id already in the table
% ------------------------------------------------------------------------------
maxId = mysql_dbquery(dbc,sprintf('SELECT MAX(%s) FROM %s',theid,theTable));
if isempty(maxId) || isempty(maxId{1}) || isnan(maxId{1})
    % No time series exist in the database yet
    maxId = 0;
else
    maxId = maxId{1}; % the maximum id -- the new items will have ids greater than this
end

% ------------------------------------------------------------------------------
%% Assemble and execute the INSERT queries
% ------------------------------------------------------------------------------
fprintf('Adding %u new %s to the %s table in %s...',...
                    sum(~isBad),theWhat,theTable,databaseName)
switch addWhat
case 'ts' % Add time series to the TimeSeries table, X at a time
          % (so as not to exceed the max_allowed_packet when transmitting the
          % time-series data).
          % An appropriate chunk size will depend on the length of time series
          % in general...
          % A kind of big problem if you add some here and then fail, because
          % then there will be inconsistencies in the Results table.
          % Safer would be to initialize results as you go in case you get a
          % fail at some chunk number... :/
    maxLength = max(TimeSeries.Length);
    % 4 at a time works for 20,000 sample time series using default mySQL settings.
    % Calibrate to this:
    numPerChunk = max([floor(20000/maxLength*4),1]); % no fewer than 1 per chunk
    numPerChunk = min([numPerChunk,500]); % no more than 500 per chunk
    SQL_AddChunked(dbc,['INSERT INTO TimeSeries (Name, Keywords, Length, ' ...
                                    'Data) VALUES'],toAdd(~isBad),numPerChunk);
case 'ops' % Add operations to the Operations table 500 at a time
    SQL_AddChunked(dbc,['INSERT INTO Operations (Name, Code, MasterLabel, ' ...
                                    'Keywords) VALUES'],toAdd(~isBad),500);
case 'mops' % Add master operations to the MasterOperations table 500 at a time
    SQL_AddChunked(dbc,['INSERT INTO MasterOperations (MasterLabel, ' ...
                                'MasterCode) VALUES'],toAdd(~isBad),500);
end
fprintf(1,' Done.\n');

% ------------------------------------------------------------------------------
% Add new entries to the Results table
% ------------------------------------------------------------------------------
if ~strcmp(addWhat,'mops')
    resultsTic = tic;
    fprintf(1,'Updating the Results table in %s...(This could take a while)...',databaseName);
    if strcmp(addWhat,'ts') && height(TimeSeries) > 1000
        fprintf(1,'\n(e.g., ~10 hours for 20,000 time series with 8,000 operations---please be patient!)...')
    end
    % The following turns out to be quite slow when adding to an existing database.
    % Maybe it does the cross join before filtering on ID. So if you have a large
    % database already, it takes ages to do the cross join, only to filter out most
    % of it anyway...?
    switch addWhat
    case 'ts'
        [~,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,op_id) SELECT t.ts_id,o.op_id ' ...
                    'FROM TimeSeries t CROSS JOIN Operations o ON t.ts_id > %u ORDER BY t.ts_id, o.op_id'],maxId));
    case 'ops'
        [~,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,op_id) SELECT t.ts_id,o.op_id ' ...
                    'FROM TimeSeries t CROSS JOIN Operations o ON o.op_id > %u ORDER BY t.ts_id, o.op_id'],maxId));
    end
    if ~isempty(emsg)
        fprintf(1,' Error. This is really really not good.\n%s',emsg);
        keyboard
    else
        fprintf(1,' initialized in %s!\n',BF_TheTime(toc(resultsTic)));
    end
end

if ~strcmp(addWhat,'mops')
    % ------------------------------------------------------------------------------
    % Update the keywords table
    % ------------------------------------------------------------------------------
    fprintf(1,'Updating the %s table in %s...',theKeywordTable,databaseName);

    % First find unique keywords from new time series by splitting against commas
    switch addWhat
    case 'ts'
        kws = TimeSeries.Keywords(~isBad);
    case 'ops'
        kws = Operations.Keywords(~isBad);
    end

    % Split into each individual keyword:
    kwSplit = regexp(kws,',','split','ignorecase');
    ukws = unique(horzcat(kwSplit{:}));
    nkw = length(ukws); % The number of unique keywords in the new set of time series
    if beVocal
        fprintf(1,'\nI found %u unique keywords in the %u new %s in %s...',...
                        nkw,sum(~isBad),theWhat,inputFile);
    end

    % How many overlap with existing keywords??:
    allkws = mysql_dbquery(dbc,sprintf('SELECT Keyword FROM %s',theKeywordTable));
    if ~isempty(allkws) % the table may be empty, in which case all keywords will be new
        isNew = ~ismember(ukws,allkws);
        % cellfun(@(x)~isempty(x),regexp(ukws,allkws,'ignorecase')); % ignore case for keywords
    else
        isNew = ones(nkw,1); % All are new
    end

    if sum(isNew) > 0
        if beVocal
            fprintf(1,['\nIt turns out that %u keywords are completely new and will be added ' ...
                        'to the %s table in %s...'],sum(isNew),theKeywordTable,databaseName);
        end
        % Add the new keywords to the Keywords table
        insertString = sprintf('INSERT INTO %s (Keyword,NumOccur) VALUES',theKeywordTable);
        fisNew = find(isNew); % Indicies of new keywords
        toAdd = cell(sum(isNew),1);
        for k = 1:sum(isNew);
            toAdd{k} = sprintf('(''%s'',0)',ukws{fisNew(k)});
        end
        SQL_AddChunked(dbc,insertString,toAdd,100);
        fprintf(1,' Added %u new keywords!\n',sum(isNew));
    else
        if beVocal
            fprintf(1,['\nIt turns out that all new keywords already exist in ' ...
                    'the %s table in %s -- there are no new keywords to add\n'],...
                    theKeywordTable,databaseName);
        end
    end

    % ------------------------------------------------------------------------------
    %% Fill new keyword relationships
    % ------------------------------------------------------------------------------
    fprintf(1,'Writing new keyword relationships to the %s table in %s...', ...
                                            theRelTable,databaseName);

    % Try doing it from scratch...:
    switch addWhat
    case 'ts'
        allNames = BF_cat(TimeSeries.Name(~isBad),',','''');
    case 'ops'
        allNames = BF_cat(Operations.Name(~isBad),',','''');
    end

    % ourids: ids matching time series or operations
    ourids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE Name IN (%s)',theid, ...
                                                        theTable,allNames));
    ourids = vertcat(ourids{:});

    % ourkids: kw_ids matching keywords
    ourkids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE Keyword IN (%s)', ...
                                            thekid,theKeywordTable,BF_cat(ukws,',','''')));
    ourkids = vertcat(ourkids{:});
    addCell = {};
    for i = 1:length(kwSplit)
        for j = 1:length(kwSplit{i})
            try
                addCell{end+1} = sprintf('(%u,%u)',ourids(i),ourkids(strcmp(kwSplit{i}{j},ukws)));
            catch
                keyboard
            end
        end
    end
    % Add them all in 1000-sized chunks
    SQL_AddChunked(dbc,sprintf('INSERT INTO %s (%s,%s) VALUES',theRelTable,theid,thekid),addCell,1000);

    % Update Nmatches in the keywords table
    fprintf(1,' Done.\nNow computing all match counts for all keywords...\n');
    SQL_FlushKeywords(addWhat);
end

% ------------------------------------------------------------------------------
%% Update links between operations and master operations
% ------------------------------------------------------------------------------
if ismember(addWhat,{'mops','ops'}) % there may be new links
    % Add mop_ids to Operations table
    fprintf(1,'Evaluating links between operations and master operations...'); tic
    updateString = ['UPDATE Operations AS o SET o.mop_id = (SELECT mop_id FROM MasterOperations AS m ' ...
                        'WHERE m.MasterLabel = o.MasterLabel) WHERE mop_id IS NULL'];
    [~,emsg] = mysql_dbexecute(dbc,updateString);

    if isempty(emsg)
        fprintf(' Done.\n');
    else
        error('\nOops! Error finding links between Operations and MasterOperations:\n%s\n',emsg);
    end

    updateString = sprintf(['UPDATE MasterOperations AS m SET NPointTo = ' ...
                    '(SELECT COUNT(o.mop_id) FROM Operations AS o WHERE m.mop_id = o.mop_id)']);
    [~,emsg] = mysql_dbexecute(dbc,updateString);
    if ~isempty(emsg)
        error('Error counting NPointTo operations for mop_id = %u\n%s\n',M_ids(k),emsg);
    end

end

% ------------------------------------------------------------------------------
%% Close database
% ------------------------------------------------------------------------------
SQL_CloseDatabase(dbc)

% ------------------------------------------------------------------------------
%% Tell the user all about it
% ------------------------------------------------------------------------------
fprintf('All tasks completed in %s.\nRead %s then added %u %s into %s.\n', ...
            BF_TheTime(toc(ticker)),inputFile,sum(~isBad),theWhat,databaseName);

if strcmp(addWhat,'ts') && ~isMatFile
    fprintf(1,['**The imported data are now in %s and no longer ' ...
                    'need to be in the Matlab path**\n'],databaseName);
end

end
