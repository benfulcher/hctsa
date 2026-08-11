function theTable = TS_ReadInputFile(addWhat,inputFile,beVocal)
% TS_ReadInputFile   Interpret a structured input file of time series, operations,
%                    or master operations into an hctsa-style table.
%
%---INPUTS:
% addWhat: 'mops' (for master operations), 'ops' (for operations), or 'ts'
%             (for time series)
% inputFile: the filename of the tab-delimited textfile to be read in [default
%            = INP_ts.txt or INP_ops.txt or INP_mops.txt]
%            The input file should be formatted with whitespace as a delimiter
%            between the entries to import.
% beVocal: if true (default) gives user feedback on the input process.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
    switch addWhat
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

if nargin < 3
    beVocal = true; % Give user feedback by default:
end

% ------------------------------------------------------------------------------
% Display welcome message:
% ------------------------------------------------------------------------------
if beVocal
    fprintf(1,'Using input file: %s.\n',inputFile);
end

% ------------------------------------------------------------------------------
% Define strings to unify the different strands of code for time series /
% operations
% ------------------------------------------------------------------------------
switch addWhat
    case 'ts'
        theWhat = 'time series';
        maxL = 50000; % the longest time series length accepted
    case 'ops'
        theWhat = 'operations';
    case 'mops'
        theWhat = 'master operations';
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

        fprintf(1,'\nHow does it look? Make sure the metadata matches the headings\n');

        % Ask the question:
        if strcmp(addWhat,'ts')
            reply = input(['If we go on, we will attempt to read all time series ' ...
                        'from file and add all ' ...
                        'data to HCTSA.mat\n<<<Type ''y'' to continue...>>> '],'s');
        else
            reply = input(['If we go on, we will attempt to add all ' ...
                        'data to HCTSA.mat\n<<<Type ''y'' to continue...>>> '],'s');
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
            drawnow
        end

        reply = input(sprintf(['Does this look ok for the first %u time series?\n' ...
                                '<<<Type ''y'' to continue...>>> '],plotNum),'s');
        if ~strcmp(reply,'y')
            fprintf(1,'I didn''t think so. Come back later...\n');
            return
        end
        close
    end

    % Ok, so we have inputData.timeSeriesData, inputData.labels, and inputData.keywords
end

% ------------------------------------------------------------------------------
% Construct a table of metadata for the time series / operations / master operations
% ------------------------------------------------------------------------------

% Prepare local IDs:
ID = (1:numItems)';

switch addWhat
case 'ts' % Prepare metadata for time series
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
                '[Note that this maximum length can be modified in TS_ReadInputFile]\n',...
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
        fprintf(1,'\nAll time-series data loaded; %u/%u passed quality tests.\n',...
                        sum(wasGood),length(wasGood));
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
    end
end

% ------------------------------------------------------------------------------
% Assemble the output table
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
fprintf(1,'Using %u %s.\n',height(theTable),theWhat);

end
