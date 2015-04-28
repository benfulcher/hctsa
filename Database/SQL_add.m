% ------------------------------------------------------------------------------
% SQL_add
% ------------------------------------------------------------------------------
% 
% Adds a set of time series, operations, or master operations to the mySQL
% database.
% 
%---INPUTS:
% addWhat: 'mops' (for master operations), 'ops' (for operations), or 'ts'
%             (for time series)
% inputFile:    the filename of the tab-delimited textfile to be read in [default
%             = INP_ts.txt or INP_ops.txt or INP_mops.txt]
%             The input file should be formatted with whitespace as a delimiter
%             between the entries to import.
%
%%---HISTORY:
% Ben Fulcher, June 2013 -- Reformulated the whole format so that only a single
% thing is uploaded at a time (ts, ops, mops), and follows a uniform and more
% transparent structure with as much overlap in syntax as possible.
%                       Added beVocal input
% Romesh Abeysuriya, Jan 2013
% Ben Fulcher, 12/1/2010: added databaseName option
% Ben Fulcher, 3/12/2009
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

function SQL_add(addWhat, inputFile, databaseName, beVocal)

% ------------------------------------------------------------------------------
%% Check inputs, set defaults:
% ------------------------------------------------------------------------------

% addWhat
% SHOULD BE 'ts', 'mops', or 'ops'
if nargin < 1 || isempty(addWhat) || ~ismember(addWhat,{'ops','ts','mops'})
    error(['Error setting first input argument -- should be ''ts'' for TimeSeries ' ...
                ', ''ops'' for Operations, or ''mops'' for Master Operations']);
end

% inputFile
if nargin < 2 || isempty(inputFile)
    % Default filenames:
    if strcmp(addWhat,'ts')
        inputFile = 'INP_ts.txt';
    else
        inputFile = 'INP_ops.txt';
    end
end
if ~exist(inputFile)
    error('Unknown file ''%s''',inputFile);
end

% databaseName
if nargin < 3
    % Use the default database (specified in sql_settings.conf) by default:
    databaseName = '';
end

% beVocal
if nargin < 4
    % Give user feedback by default:
    beVocal = 1;
end

% ------------------------------------------------------------------------------
% Display welcome message:
% ------------------------------------------------------------------------------
if beVocal,
    fprintf(1,'Using input file: %s\n',inputFile);
else
    fprintf(1,'No user feedback shown (use beVocal to get information about your import).\n');
end
ticker = tic;

% ------------------------------------------------------------------------------
%% Open Database
% ------------------------------------------------------------------------------
[dbc, databaseName] = SQL_opendatabase(databaseName);

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
        thektable = 'TimeSeriesKeywords';
        theRelTable = 'TsKeywordsRelate';
        thename = 'Filename';
        maxL = 50000; % the longest time series length accepted in the database
    case 'ops'
        theWhat = 'operations';
        theid = 'op_id';
        thekid = 'opkw_id';
        theTable = 'Operations';
        thektable = 'OperationKeywords';
        theRelTable = 'OpKeywordsRelate';
        thename = 'OpName';
    case 'mops'
        theWhat = 'master operations';
        theid = 'mop_id';
        theTable = 'MasterOperations';
end

% ------------------------------------------------------------------------------
% Check that the table exists in the datbase
% ------------------------------------------------------------------------------
existString = ['SHOW TABLES LIKE ''' theTable ''''];
[output,emsg] = mysql_dbquery(dbc,existString);
if isempty(output) % Table doesn't exist
    error(['Table %s doesn''t exist in the database %s.\n' ...
                'Use: (1) install.m to set the database system up from scratch,\n' ...
                '(2) SQL_create_all_tables to create empty tables that can later be filled with custom libraries of operations,' ...
                '(3) SQL_reset to drop all tables in the database, and repopulate all tables with the default library of operations.'],...
                databaseName);
end

% ------------------------------------------------------------------------------
%% Open and read the input file
% ------------------------------------------------------------------------------

% Determine if it's a .mat file:
if strcmp(inputFile(end-3:end),'.mat');
    isMatFile = 1;
else
    isMatFile = 0;
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
            fprintf(1,['Need to format %s (Time Series input file) as: Filename ' ...
                                                                'Keywords\n'],inputFile)
            fprintf(1,'Assuming no header line\n')
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
        end
    	dataIn = textscan(fid,'%s %s','CommentStyle','#','CollectOutput',1);
    
    case 'ops' % Read the operations input file:
        if beVocal
            fprintf(1,['Need to format %s (Operations input file) as: OperationCode ' ...
                                            'OperationName OperationKeywords\n'],inputFile)
            fprintf(1,'Assuming no header lines\n')
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
        end
        dataIn = textscan(fid,'%s %s %s','CommentStyle','#','CollectOutput',1);    
    
    case 'mops' % Read the master operations input file:
        if beVocal
            fprintf(1,'Need to format %s (Master Operations input file) as: MasterCode MasterLabel\n',inputFile)
            fprintf(1,'Assuming no header lines\n')
            fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
            fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
        end
        dataIn = textscan(fid,'%s %s','CommentStyle','#','CollectOutput',1);
    end

    fclose(fid);
    

    % ------------------------------------------------------------------------------
    % Show the user what's been imported:
    % ------------------------------------------------------------------------------
    dataIn = dataIn{1}; % Collect one big matrix of cells
    numItems = size(dataIn,1); % Number of items in the input file

    if numItems == 0, error('The input file ''%s'' seems to be empty??',inputFile), end

    if beVocal
        fprintf(1,'Found %u %s in %s, I think. Take a look:\n\n',numItems,theWhat,inputFile)
        switch addWhat
        case 'ts'
            fprintf(1,'%s\t%s\n','-Filename-','-Keywords-')
            fprint_ts = @(x) fprintf('%s\t%s\n',dataIn{x,1},dataIn{x,2});
        case 'ops'
            fprintf(1,'%s\t%s\t%s\n','-Operation Name-','-Operation Code-','-Operation Keywords-')
            fprint_ops = @(x) fprintf('%s\t%s\t%s\n',dataIn{x,1},dataIn{x,2},dataIn{x,3});
        case 'mops'
            fprintf(1,'%s\t%s\n','-Master Code-','-Master Label-')
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
            fprintf(1,'..................(%u).....................\n',max(numItems-6,0))
            for i = max(numItems-2,4):numItems
                switch addWhat
                case 'ts', fprint_ts(i);
                case 'ops', fprint_ops(i);
                case 'mops', fprint_mops(i);
                end
            end
        end
    
        fprintf(1,['\nHow does it look? Make sure the time series and everything ' ...
                                                'match up to their headings\n'])

        reply = input(['If we go on, we will attempt to read all timeseries ' ...
                        'from file and add all ' ...
                        'data to the database.\n<<<Type ''y'' to continue...>>>'],'s');
    
        if ~strcmp(reply,'y')
            fprintf(1,'I didn''t think so. Come back later...\n')
            return
        end
    end

    fprintf(1,'%s read.\n',inputFile)
    
else
    % Only allow .mat file input when importing time series
    if ~strcmp(addWhat,'ts')
        error(['.mat file input type only supported for importing time series. ' ...
                'Please specify a text file for importing master operations or operations.']);
    end
        
    % Load the 3 cells specifying the data:
    if ~exist(inputFile)
        error('Could not load specified input .mat file: %s',inputFile);
    end
    inputData = load(inputFile,'timeSeriesData','labels','keywords');
    
    % Check that they're all as they should be:
    if ~isfield(inputData,'timeSeriesData') || ~isfield(inputData,'labels') || ~isfield(inputData,'keywords') ...
                || ~(iscell(inputData.timeSeriesData) || isnumeric(inputData.timeSeriesData)) ...
                || ~iscell(inputData.labels) || ~iscell(inputData.keywords)
        error(['Expecting input file, %s, to contain: ''timeSeriesData'' (cell or matrix), ' ...
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
                                    '(loaded from %s).'],numItems,inputFile)
        fprintf(1,['Will store time-series data from matlab file in the database to an' ...
                        ' accuracy of 6 significant figures...\n']);
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
        
        reply = input(sprintf(['Does this look ok for the first %u time series?\n If we continue, ' ...
                        'we will attempt to add all %u time series in the input file' ...
                        ' to the database.\n<<<Type ''y'' to continue...>>>'],plotNum,numItems),'s');
        if ~strcmp(reply,'y')
            fprintf(1,'I didn''t think so. Come back later...\n')
            return
        end
        close
    end
    
    % Ok, so we have inputData.timeSeriesData, inputData.labels, and inputData.keywords
end

esc = @RA_sqlescapestring; % Inline function to add escape strings to format mySQL queries

% ------------------------------------------------------------------------------
% Construct a more intuitive structure array for the time series / operations /
% master operations Fill a cell, toAdd, containing mySQL INSERT commands for
% each item in the input file.
% ------------------------------------------------------------------------------
if beVocal
    fprintf(1,['Preparing mySQL statements to add %u %s to the ' ...
                                'database %s...'],numItems,theWhat,databaseName);
end
toAdd = cell(numItems,1);
switch addWhat
case 'ts' % Prepare toAdd cell for time series
    if beVocal; figure('color','w','WindowStyle','docked'); end
    for j = 1:numItems
        
        % Assign filename and keywords strings to this time series, and load it as x
        if isMatFile
            timeseries(j).Filename = inputData.labels{j};
            timeseries(j).Keywords = inputData.keywords{j}; % Take out inverted commas from keywords lists
            if iscell(inputData.timeSeriesData)
                x = inputData.timeSeriesData{j};
            else
                x = inputData.timeSeriesData(j,:);
            end
        else
            timeseries(j).Filename = dataIn{j,1};
            timeseries(j).Keywords = regexprep(dataIn{j,2},'\"',''); % Take out inverted commas from keywords lists
            % Read the time series from its filename:
            try
                x = dlmread(timeseries(j).Filename);
            catch emsg
                fprintf(1,'%s\n',emsg.message)
                error(['\nCould not read the data file for ''%s''.' ...
                                'Check that it''s in Matlab''s path.'], ...
                                    timeseries(j).Filename)
            end
        end

        % Assign the length of the time series
        timeseries(j).Length = length(x);
        
        % Check if the time series contains any NaN of Inf values:
        if any(isnan(x)) || any(~isfinite(x))
            fprintf(1,['\nDid you know that the time series %s contains special values' ...
                        ' (e.g., NaN or Inf)...?\n'],which(timeseries(j).Filename))
            fprintf(1,'I''m not quite sure what to do with this... Please reformat.\n')
            return
        end
        
        % If this time series is longer than the maximum allowed, then exit:
        if length(x) > maxL
            fprintf(['\n%s contains %u samples, this framework can efficiently ' ...
                            'deal with time series up to %u samples\n'],...
                                timeseries(j).Filename,timeseries(j).Length,maxL)
            fprintf(1,['Safest not to do anything now -- perhaps you can remove this ' ...
                                    'time series from the input file?']); return
        end
        
        % Now assign the time-series data (stored as singles in the case of a .mat file data)
        if isMatFile
            xtext = sprintf('%.6g,',x); % keeps 6 figures of accuracy
            xtext = xtext(1:end-1); % remove trailing comma
        else
            % Read in the time series from file as strings using textscan
            fid = fopen(timeseries(j).Filename);
            timeSeriesData_text = textscan(fid,'%s');
            fclose(fid);
        
            % Turn the time series into a comma-delimited string
            timeSeriesData_text = timeSeriesData_text{1}; % cell of time series values as strings
            xtext = ''; % make a comma-delimited text version of the time series as xtext
            for k = 1:length(timeSeriesData_text)
                xtext = [xtext,',',timeSeriesData_text{k}];
            end
            xtext = xtext(2:end);
        end
        timeseries(j).Data = xtext;
        
        % ------------------------------------------------------------------------------
        % Prepare the data to be added to the database in an INSERT command:
        toAdd{j} = sprintf('(''%s'',''%s'',%u,''%s'')', ...
                            esc(timeseries(j).Filename),esc(timeseries(j).Keywords), ...
                            timeseries(j).Length,timeseries(j).Data);
        
        if beVocal % plot the time series
            nsubplots = min(numItems,5);
            subplot(nsubplots,1,mod(j-1,nsubplots)+1);
            plot(x,'-k'); xlim([1,length(x)]);
            titletext = sprintf('\n[%u/%u] %s (%u), keywords = %s',j,numItems,timeseries(j).Filename,timeseries(j).Length,timeseries(j).Keywords);
            title(titletext,'interpreter','none');
            fprintf(1,'%s --- loaded successfully.',titletext)
            pause(0.01); % wait 0.01 second to show the plotted time series
        end
    end
    if beVocal, fprintf(1,'\nAll time-series data loaded, ready to be uploaded to the mySQL database.\n'); end

    % Check for duplicates in the input file:
    if length(unique({timeseries.Filename})) < length(timeseries)
        error('Input file contains duplicates.');
    end
    
case 'mops' % Prepare toAdd cell for master operations
    for j = 1:numItems
        master(j).MasterCode = dataIn{j,1};
        master(j).MasterLabel = dataIn{j,2};
        toAdd{j} = sprintf('(''%s'', ''%s'')',esc(master(j).MasterLabel),esc(master(j).MasterCode));
    end
    if beVocal, fprintf(1,' Done.\n'); end
    
case 'ops' % Prepare toAdd cell for operations        
    for j = 1:numItems
        operation(j).Code = dataIn{j,1};
        operation(j).Name = dataIn{j,2};
        operation(j).Keywords = dataIn{j,3};
        if strfind(operation(j).Code,'(') % single operation
            operation(j).MasterLabel = '';
            toAdd{j} = sprintf('(''%s'', ''%s'',NULL,''%s'')',esc(operation(j).Name),esc(operation(j).Code),esc(operation(j).Keywords));
        else % pointer operation
            operation(j).MasterLabel = strtok(operation(j).Code,'.');
            toAdd{j} = sprintf('(''%s'', ''%s'',''%s'',''%s'')',esc(operation(j).Name),esc(operation(j).Code),esc(operation(j).MasterLabel),esc(operation(j).Keywords));
        end
    end
    if beVocal, fprintf(1,' Done.\n'); end
    
    % Check for duplicates in the input file:
    [uniqueOpNames,ia] = unique({operation.Name});
    if length(uniqueOpNames) < length(operation)
        warning(['Input file contains %u duplicate entries, which are being removed.\n' ...
                        'Inputting %u -> %u operations...'], ...
            length(operation)-length(uniqueOpNames),length(operation),length(uniqueOpNames));
        % Only keep the unique ones:
        operation = operation(ia);
        toAdd = toAdd(ia);
        numItems = length(operation);
        fprintf(1,'We now have %u operations to input...\n',numItems);
    end    
end

% ------------------------------------------------------------------------------
%% Check for duplicates
% ------------------------------------------------------------------------------
if beVocal, fprintf(1,'Checking for duplicates already in the database... '); end
switch addWhat
case 'ts'
    existingFilenames = mysql_dbquery(dbc,sprintf('SELECT FileName FROM TimeSeries'));
    isDuplicate = ismember({timeseries.Filename},existingFilenames); % isDuplicate = 1 if the item already exists
case 'ops'
    existingOperationNames = mysql_dbquery(dbc,'SELECT OpName FROM Operations');
    isDuplicate = ismember({operation.Name},existingOperationNames); % isDuplicate = 1 if the MasterLabel already exists    
case 'mops'
    existing = mysql_dbquery(dbc,'SELECT MasterLabel FROM MasterOperations');
    isDuplicate = ismember({master.MasterLabel},existing); % isDuplicate = 1 if the MasterLabel already exists
end
if beVocal, fprintf(1,'done.\n'); end

% Tell the user about duplicates
if all(isDuplicate)
    fprintf(1,'All %u %s from %s already exist in %s---no new %s to add!\n',numItems,theWhat,inputFile,databaseName,theWhat);
    return
elseif sum(isDuplicate) > 0
    if beVocal
        fprintf(1,'I found %u duplicate %s already in the database %s!\n',sum(isDuplicate),theWhat,databaseName)
        fprintf(1,'There are %u new %s to add to %s.\n',sum(~isDuplicate),theWhat,databaseName)
    end
end

% ------------------------------------------------------------------------------
%% Select the maximum id already in the table
% ------------------------------------------------------------------------------
maxid = mysql_dbquery(dbc,sprintf('SELECT MAX(%s) FROM %s',theid,theTable));
if isempty(maxid) || isempty(maxid{1}) || isnan(maxid{1})
    % No time series exist in the database yet
    maxid = 0;
else
    maxid = maxid{1}; % the maximum id -- the new items will have ids greater than this
end

% ------------------------------------------------------------------------------
%% Assemble and execute the INSERT queries
% ------------------------------------------------------------------------------
fprintf('Adding %u new %s to the %s table in %s...',sum(~isDuplicate),theWhat,theTable,databaseName)
switch addWhat
case 'ts' % Add time series to the TimeSeries table just 5 at a time
          % (so as not to exceed the max_allowed_packet when transmitting the 
          % time-series data) Appropriate chunk size will depend on the length of
          % time series in general

    SQL_add_chunked(dbc,['INSERT INTO TimeSeries (FileName, Keywords, Length, ' ...
                                    'Data) VALUES'],toAdd,isDuplicate,5);
case 'ops' % Add operations to the Operations table 500 at a time
    SQL_add_chunked(dbc,['INSERT INTO Operations (OpName, Code, MasterLabel, ' ...
                                    'Keywords) VALUES'],toAdd,isDuplicate,500);        
case 'mops' % Add master operations to the MasterOperations table 500 at a time
    SQL_add_chunked(dbc,['INSERT INTO MasterOperations (MasterLabel, ' ...
                                'MasterCode) VALUES'],toAdd,isDuplicate,500);
end
fprintf(1,' done.\n')

% ------------------------------------------------------------------------------
% Add new entries to the Results table
% ------------------------------------------------------------------------------
if ~strcmp(addWhat,'mops')
    resultstic = tic;
    if beVocal
        fprintf(1,'Updating the Results Table in %s (this could take a while, please be patient!)...',databaseName)
    end
    switch addWhat
    case 'ts'
        [~,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,op_id) SELECT t.ts_id,o.op_id ' ...
                    'FROM TimeSeries t CROSS JOIN Operations o ON t.ts_id > %u ORDER BY t.ts_id, o.op_id'],maxid));
    case 'ops'
        [~,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,op_id) SELECT t.ts_id,o.op_id ' ...
                    'FROM TimeSeries t CROSS JOIN Operations o ON o.op_id > %u ORDER BY t.ts_id, o.op_id'],maxid));
    end
    if ~isempty(emsg),
        fprintf(1,' error. This is really really not good.\n');
        keyboard
    else
        if beVocal, fprintf(1,' initialized in %s!\n',BF_thetime(toc(resultstic))); end
    end
end

if ~strcmp(addWhat,'mops')
    % ------------------------------------------------------------------------------
    % Update the keywords table
    % ------------------------------------------------------------------------------
    fprintf(1,'Updating the %s table in %s...',thektable,databaseName)

    % First find unique keywords from new time series by splitting against commas
    switch addWhat
    case 'ts'
        kws = {timeseries(~isDuplicate).Keywords};
    case 'ops'
        kws = {operation(~isDuplicate).Keywords};
    end
    
    kwsplit = cell(length(kws),1); % Split into each individual keyword
    ukws = {};
    for i = 1:length(kws)
        kwsplit{i} = regexp(kws{i},',','split','ignorecase');
        for j = 1:length(kwsplit{i})
            if ~ismember(kwsplit{i}{j},ukws) % add it to ukws
                ukws{end+1} = kwsplit{i}{j};
            end
        end
    end
    nkw = length(ukws); % The number of unique keywords in the new set of time series
    if beVocal, fprintf(1,'\nI found %u unique keywords in the %u new %s in %s...',nkw,sum(~isDuplicate),theWhat,inputFile); end

    % How many overlap with existing keywords??:
    allkws = mysql_dbquery(dbc,sprintf('SELECT Keyword FROM %s',thektable));
    if ~isempty(allkws) % the table may be empty, in which case all keywords will be new
        isnew = ~ismember(ukws,allkws);
        % cellfun(@(x)~isempty(x),regexp(ukws,allkws,'ignorecase')); % ignore case for keywords
    else
        isnew = ones(nkw,1); % All are new
    end
    
    if sum(isnew) > 0
        if beVocal
            fprintf(1,['\nIt turns out that %u keywords are completely new and will be added ' ...
                        'to the %s table in %s...'],sum(isnew),thektable,databaseName)
        end
        % Add the new keywords to the Keywords table
        insertstring = sprintf('INSERT INTO %s (Keyword,NumOccur) VALUES',thektable);
        toAdd = cell(sum(isnew),1);
        fisnew = find(isnew); % Indicies of new keywords
        for k = 1:sum(isnew);
            toAdd{k} = sprintf('(''%s'',0)',ukws{fisnew(k)});
        end
        SQL_add_chunked(dbc,insertstring,toAdd);
        fprintf(1,' added %u new keywords!\n',sum(isnew))
    else
        if beVocal
            fprintf(1,['\nIt turns out that all new keywords already exist in ' ...
                        'the %s table in %s -- there are no new keywords to add\n'],sum(isnew),thektable,databaseName)
        end
    end
    
    % ------------------------------------------------------------------------------
    %% Fill new keyword relationships
    % ------------------------------------------------------------------------------
    fprintf(1,'Writing new keyword relationships to the %s table in %s...', ...
                                            theRelTable,databaseName)

    % Try doing it from scratch...:
    switch addWhat
    case 'ts'
        allNames = BF_cat({timeseries(~isDuplicate).Filename},',','''');
    case 'ops'
        allNames = BF_cat({operation(~isDuplicate).Name},',','''');
    end
    ourids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE %s IN (%s)',theid, ...
                                                        theTable,thename,allNames));
    ourids = vertcat(ourids{:}); % ids matching FileNames/OpNames
    ourkids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE Keyword IN (%s)', ...
                                            thekid,thektable,BF_cat(ukws,',','''')));
    ourkids = vertcat(ourkids{:}); % ids matching FileNames/OpNames
    % nkwrels = sum(cellfun(@(x)length(x),kwsplit)); % number of keyword relationships in the input file
    addCell = {};
    for i = 1:length(kwsplit)
        for j = 1:length(kwsplit{i})
            addCell{end+1} = sprintf('(%u,%u)',ourids(i),ourkids(strcmp(kwsplit{i}{j},ukws)));
        end
    end
    SQL_add_chunked(dbc,sprintf('INSERT INTO %s (%s,%s) VALUES',theRelTable,theid,thekid),addCell); % add them all in chunks
        
    % Increment Nmatches in the keywords table
    fprintf(1,' done.\nNow calculating the match counts for keywords...')
    % Redo them from scratch should be easier actually...?
    for k = 1:nkw % keywords implicated in this import
        SelectString = sprintf('(SELECT %s FROM %s WHERE Keyword = ''%s'')',thekid,thektable,ukws{k});
        theopkw = mysql_dbquery(dbc,SelectString);
        updateString = sprintf('UPDATE %s SET NumOccur = (SELECT COUNT(*) FROM %s WHERE %s = %u) WHERE %s = %u', ...
                                    thektable,theRelTable,thekid,theopkw{1},thekid,theopkw{1});
        [~,emsg] = mysql_dbexecute(dbc, updateString);
        if ~isempty(emsg)
            error('\n Error updating keyword count in %s\n%s',thektable,emsg)
        end
    end
    % for k = 1:nkw % for each unique keyword in the keyword table...
    %     % nnkw = sum(cellfun(@(x)ismember(ukws{k},x),kwsplit));
    %     Selectopkwid = sprintf('(SELECT %s FROM %s WHERE Keyword = ''%s'')',thekid,thektable,ukws{k});
    %     SelectCount = sprintf(['SELECT COUNT(*) FROM %s WHERE %s = %s ' ...
    %                             'AND %s > %u'],theRelTable,thekid,Selectopkwid,theid,maxid);
    %     updateString = sprintf(['UPDATE %s SET NumOccur = NumOccur + (%s) ' ...
    %                             'WHERE Keyword = ''%s'''],thektable,SelectCount,ukws{k});
    %     [~,emsg] = mysql_dbexecute(dbc, updateString);
    %     if ~isempty(emsg)
    %         keyboard
    %         fprintf(1,'\n Error updating keyword count in %s',thektable)
    %     end
    % end
    fprintf(1,' done.\n')
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
    
    % ---------------------
    % Should probably add a check that the update happened; i.e., that no mop_id are left NULL
    % after this process (e.g., if no Master Operations exist, then this will complete with no
    % error; even though no mop_ids were actually assigned...)

    
    %----ALTERNATIVE: FILL A LINKING TABLE AS A JOIN ON MasterLabel:
    % Delete the linking table and recreate it from scratch is easiest
    % fprintf(1,'Filling MasterPointerRelate...');
    % mysql_dbexecute(dbc,'DELETE FROM MasterPointerRelate');
    % InsertString = ['INSERT INTO MasterPointerRelate SELECT m.mop_id,o.op_id FROM MasterOperations ' ...
    %                         'm JOIN Operations o ON m.MasterLabel = o.MasterLabel'];
    % [~,emsg] = mysql_dbexecute(dbc,InsertString);

    if isempty(emsg)
        fprintf(' done.\n');
    else
        error('\nOops! Error finding links between Operations and MasterOperations:\n%s\n',emsg);
    end
    
    %     % if strcmp(addWhat,'ops')
    %     %     % operations were imported -- match their MasterLabels with elements of the MasterOperations table using mySQL JOIN
    %     %     InsertString = ['INSERT INTO MasterPointerRelate SELECT m.mop_id,o.op_id FROM MasterOperations m JOIN ' ...
    %     %                         'Operations o ON m.MasterLabel = o.MasterLabel WHERE o.op_id > %u',maxid];
    %     % else
    %     %     InsertString = ['INSERT INTO MasterPointerRelate SELECT m.mop_id,o.op_id FROM MasterOperations m JOIN ' ...
    %     %                         'Operations o ON m.MasterLabel = o.MasterLabel WHERE m.mop_id > %u',maxid];
    %     % end
    %     
    
    updateString = sprintf(['UPDATE MasterOperations AS m SET NPointTo = ' ...
                    '(SELECT COUNT(o.mop_id) FROM Operations AS o WHERE m.mop_id = o.mop_id)']);
    [~,emsg] = mysql_dbexecute(dbc, updateString);
    if ~isempty(emsg)
        error('Error counting NPointTo operations for mop_id = %u\n%s\n',M_ids(k),emsg);
    end
    
    % M_ids = mysql_dbquery(dbc,'SELECT mop_id FROM MasterOperations');
    % M_ids = vertcat(M_ids{:}); % vector of master_ids    
    % for k = 1:length(M_ids)
    %     updateString = sprintf(['UPDATE MasterOperations SET NPointTo = ' ...
    %                     '(SELECT COUNT(mop_id) FROM Operations WHERE mop_id = %u)' ...
    %                         'WHERE mop_id = %u'],M_ids(k),M_ids(k));
    %     [~,emsg] = mysql_dbexecute(dbc, updateString);
    %     if ~isempty(emsg)
    %         fprintf(1,'Error counting NPointTo operations for mop_id = %u\n',M_ids(k));
    %         fprintf(1,'%s\n',emsg)
    %         keyboard
    %     end
    % end

end

% ------------------------------------------------------------------------------
%% Close database
% ------------------------------------------------------------------------------
SQL_closedatabase(dbc)

% ------------------------------------------------------------------------------
%% Tell the user all about it
% ------------------------------------------------------------------------------
fprintf('All tasks completed in %s.\nRead %s then added %u %s into %s.\n', ...
            BF_thetime(toc(ticker)), inputFile, sum(~isDuplicate), theWhat, databaseName);

if strcmp(addWhat,'ts')
    fprintf(1,['**The imported data are now in %s and no longer ' ...
                    'need to be in the Matlab path.\n'],databaseName);
end

end