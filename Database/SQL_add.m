function SQL_add(importwhat, INPfile, dbname, bevocal)
% Adds a time series or operation to the database with the specified columns, collabs, of the table given
% INPUTS: importwhat -- 'mops', 'ops', or 'ts' -- write to which of these?
% 		  inpfilename -- the filename of the tab-delimited textfile to be read in [default = INP_ts.txt or INP_ops.txt]
% The input file should be formatted with tabs (\t) as delimiters between the entries, and in the order specified below
% the column labels of the TimeSeries or Operations table in the database are in this order:
% 			 'FileName', 'Keywords', 'Quantity', 'Unit', 'SamplingRate' for TimeSeries
% 			 'Pointer'(M,S,P),'Code','OpName','Keywords' for Operations (with keywords only for children)
% Ben Fulcher 3/12/2009
% Ben Fulcher 12/1/2010: added dbname option
% Romesh Jan 2013
% Ben Fulcher June 2013 -- reformulated the whole format so that only a single thing is uploaded at a time (ts, ops, mops), and follows a uniform and more transparent structure with as much overlap in syntax as possible. Added bevocal input

%% CHECK INPUTS:
% SHOULD BE TS, MOP, or OP -- or can iterate through each possibility
% Nice to make code that inports a given type of thing into a given table
% importwhat = 'ts'; INPfile = 'INP_ts.txt'; dbname = '';

% importwhat
if nargin < 1 || isempty(importwhat) || ~ismember(importwhat,{'ops','ts','mops'})
    error('Error setting first input argument -- should be ''ts'' for TimeSeries or ''ops'' for Operations');
end

% inpfilename
if nargin < 2 || isempty(INPfile)
    if strcmp(importwhat,'ts')
        INPfile = 'INP_ts.txt';
    else
        INPfile = 'INP_ops.txt';
    end
end

% dbname
if nargin < 3
    dbname = ''; % use the default database specified in SQL_opendatabase
end

% bevocal
if nargin < 4
    bevocal = 1; % gives lots and lots of user feedback
end

if bevocal, fprintf(1,'Using input file %s',INPfile); end
ticker = tic;

%% Open Database
[dbc, dbname] = SQL_opendatabase(dbname);

% Define strings to unify the different strands of code for time series / operations
switch importwhat
    case 'ts'
        thewhat = 'time series';
        theid = 'ts_id';
        thekid = 'tskw_id';
        thetable = 'TimeSeries';
        thektable = 'TimeSeriesKeywords';
        thereltable = 'TsKeywordsRelate';
        thename = 'Filename';
        maxL = 40000; % the longest time series length for the database
    case 'ops'
        thewhat = 'operations';
        theid = 'm_id';
        thekid = 'mkw_id';
        thetable = 'Operations';
        thektable = 'OperationKeywords';
        thereltable = 'OpKeywordsRelate';
        thename = 'OpName';
    case 'mops'
        thewhat = 'master operations';
        theid = 'mop_id';
        thetable = 'MasterOperations';
        
end


%% 1. Open, read the input file
fid = fopen(INPfile);
switch importwhat
case 'ts' % Read the time series input file:
    if bevocal
        fprintf(1,'Need to format %s (Time Series input file) as: Filename Keywords\n',INPfile)
        fprintf(1,'Assuming no header line\n')
        fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
        fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
    end
	datain = textscan(fid,'%s %s','CommentStyle','%','CollectOutput',1); % 'HeaderLines',1,
case 'ops' % Read the operations input file:
    if bevocal
        fprintf(1,'Need to format %s (Operations input file) as: OperationName OperationCode OperationKeywords\n',INPfile)
        fprintf(1,'Assuming no header lines\n')
        fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
        fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
    end
    datain = textscan(fid,'%s %s %s','CommentStyle','%','CollectOutput',1);    
case 'mops' % Read the master operations input file:
    if bevocal
        fprintf(1,'Need to format %s (Master Operations input file) as: MasterCode MasterLabel\n',INPfile)
        fprintf(1,'Assuming no header lines\n')
        fprintf(1,'Use whitespace as a delimiter and \\n for new lines...\n')
        fprintf(1,'(Be careful that no additional whitespace is in any fields...)\n')
    end
    datain = textscan(fid,'%s %s','CommentStyle','%','CollectOutput',1);
end
fclose(fid);


% Show the user what's been imported:
datain = datain{1}; % collect one big matrix of cells
nits = size(datain,1); % number of items in the input file

if nits == 0, error(['Input file ' INPfile ' seems to be empty??']), end

if bevocal
    fprintf(1,'Found %u %s in %s, I think. Take a look:\n',nits,thewhat,INPfile)
    switch importwhat
    case 'ts'
        fprintf(1,'%s\t%s\n','Filename','Keywords')
        fprint_ts = @(x) fprintf('%s\t%s\n',datain{x,1},datain{x,2});
    case 'ops'
        fprintf(1,'%s\t%s\t%s\n','Operation Name','Operation Code','Operation Keywords')
        fprint_ops = @(x) fprintf('%s\t%s\t%s\n',datain{x,1},datain{x,2},datain{x,3});
    case 'mops'
        fprintf(1,'%s\t%s\n','Master Code','Master Label')
        fprint_mops = @(x) fprintf('%s\t%s\n',datain{x,1},datain{x,2});
    end
    for i = 1:min(3,nits)
        switch importwhat
        case 'ts', fprint_ts(i);
        case 'ops', fprint_ops(i);
        case 'mops', fprint_mops(i);
        end
    end
    if nits > 3
        fprintf(1,'..................(%u).....................\n',nits-3)
        for i = max(nits-2,4):nits
            switch importwhat
            case 'ts', fprint_ts(i);
            case 'ops', fprint_ops(i);
            case 'mops', fprint_mops(i);
            end
        end
    end
    fprintf(1,'How does it look? Make sure the time series and everything match up to their headings\n')
    reply = input('If we go on, we will attempt to read all timeseries from file and add all data to the database. Continue...? [y]','s');
    if ~strcmp(reply,'y')
        fprintf(1,'I didn''t think so. Come back later...\n')
        return
    end
end
fprintf(1,'%s read.\n',INPfile)
esc = @sqlescapestring; % inline function to add escape strings to format mySQL queries

% Construct a more intuitive structure array for the time series / operations / master operations
% Fill a cell, toadd, containing mySQL INSERT commands for each item in the input file:
if bevocal
    fprintf(1,'Preparing mySQL INSERT statements to add %u %s to the database %s\n',nits,thewhat,dbname);
end
toadd = cell(nits,1);
resave = 0; % need user permission to save over existing time series
switch importwhat
case 'ts' % Prepare toadd cell for time series
    if bevocal; figure('color','w','WindowStyle','docked'); end
        for j = 1:nits
        timeseries(j).Filename = datain{j,1};
        timeseries(j).Keywords = regexprep(datain{j,2},'\"',''); % Take out inverted commas from keywords lists

        % Read the time series to record its length:
        try
            x = dlmread(timeseries(j).Filename);
            timeseries(j).Length = length(x);
            % if size(x,2) > size(x,1); % a row vector
            %     if size(x,1)==1
            %         % Get permission once
            %         if resave == 0
            %             fprintf(1,['Looks like some time series files (%s) are row vectors instead of column vectors.' ...
            %                             'Can I resave over them?\n'],which(timeseries(j).Filename))
            %             reply = input('''y'' for yes (will remember this answer for other time series), or any other key to quit','s');
            %             if strcmp(reply,'y')
            %                 resave = 1;
            %             else
            %                 fprintf(1,'Exiting. You''ll have to reformat some of the time series in %s.\n',INPfile)
            %                 return
            %             end
            %         end
            %         x = x';
            %         dlmwrite(which(timeseries(j).Filename),x);
            %         fprintf(1,['%s is a row vector -- has been transposed and resaved as ''%s'' ' ...
            %                         'as a column vector\n'],timeseries(j).Filename,which(timeseries(j).Filename))
            %     else
            %         fprintf(1,['%s is not a row or column vector -- not sure how to ' ...
            %                     'read this as a time series...!\n'],which(timeseries(j).Filename))
            %         fprintf(1,'Maybe check and reformat the time series in %s...? Exiting.\n',INPfile)
            %         return
            %     end
            % end
            
            if any(isnan(x)) || any(~isfinite(x))
                fprintf(1,['Did you know that the time series %s contains special values' ...
                            ' (e.g., NaN or Inf)...?\n'],which(timeseries(j).Filename))
                fprintf(1,'I''m not quite sure what to do with this... Please reformat.\n')
                return
            end
            
            if length(x) > maxL % this time series is too long -- exit
                fprintf(['%s contains %u samples, this framework can efficiently ' ...
                                'deal with time series up to %u samples\n'],timeseries(j).Filename,timeseries(j).Length,maxL)
                fprintf(1,'Safest not to do anything now -- perhaps you can remove this time series from the input file?'); return
            end
            % Read in the time series from file as strings using textscan
            fid = fopen(timeseries(j).Filename); tstext = textscan(fid,'%s'); fclose(fid);
            % Turn into one long comma-delimited string
            tstext = tstext{1}; % cell of time series values as strings
            xtext = ''; % make a comma-delimited text version of the time series as xtext
            for k = 1:length(tstext)
                xtext = [xtext,',',tstext{k}];
            end
            xtext = xtext(2:end);
            timeseries(j).Data = xtext;

        catch emsg
            fprintf(1,'%s\n',emsg)
            error(sprintf(['Could not read the data file for ''%s''.' ...
                                    'Check that it''s in Matlab''s path.'],timeseries(j).Filename))
        end
        toadd{j} = sprintf('(''%s'',''%s'',%u,''%s'')',esc(timeseries(j).Filename),esc(timeseries(j).Keywords),timeseries(j).Length,timeseries(j).Data);
        
        if bevocal % plot the time series
            nsubplots = min(nits,5);
            subplot(nsubplots,1,mod(j-1,nsubplots)+1);
            plot(x,'-k'); xlim([1,length(x)]);
            titletext = sprintf('[%u/%u] %s (%u), keywords = %s --- read',j,nits,timeseries(j).Filename,timeseries(j).Length,timeseries(j).Keywords);
            title(titletext,'interpreter','none');
            fprintf(1,titletext)
            pause(0.2); % wait 0.2 seconds
        end
    end
    
case 'mops' % Prepare toadd cell for master operations
    for j = 1:nits
        master(j).MasterCode = datain{j,1};
        master(j).MasterLabel = datain{j,2};
        toadd{j} = sprintf('(''%s'', ''%s'')',esc(master(j).MasterLabel),esc(master(j).MasterCode));
    end
    
case 'ops' % Prepare toadd cell for operations        
    for j = 1:nits
        operation(j).Code = datain{j,1};
        operation(j).Name = datain{j,2};
        operation(j).Keywords = datain{j,3};
        if strfind(operation(j).Code,'(') % single operation
            operation(j).MasterLabel = '';
            toadd{j} = sprintf('(''%s'', ''%s'',NULL,''%s'')',esc(operation(j).Name),esc(operation(j).Code),esc(operation(j).Keywords));
        else % pointer operation
            operation(j).MasterLabel = strtok(operation(j).Code,'.');
            toadd{j} = sprintf('(''%s'', ''%s'',''%s'',''%s'')',esc(operation(j).Name),esc(operation(j).Code),esc(operation(j).MasterLabel),esc(operation(j).Keywords));
        end
        
    end
end
if bevocal, fprintf(1,'done.\n'); end


% Check for duplicates
if bevocal, fprintf(1,'Checking for duplicates already in the database... '); end
switch importwhat
case 'ts'
    ExistingFilenames = mysql_dbquery(dbc,sprintf('SELECT FileName FROM TimeSeries'));
    isduplicate = ismember({timeseries.Filename},ExistingFilenames); % isduplicate = 1 if the item already exists
case 'ops'
    ExistingOperationNames = mysql_dbquery(dbc,'SELECT OpName FROM Operations');
    isduplicate = ismember({operation.Name},ExistingOperationNames); % isduplicate = 1 if the MasterLabel already exists    
case 'mops'
    existing = mysql_dbquery(dbc,'SELECT MasterLabel FROM MasterOperations');
    isduplicate = ismember({master.MasterLabel},existing); % isduplicate = 1 if the MasterLabel already exists
end
if bevocal, fprintf(1,'done.\n'); end

% Tell the user about duplicates

if all(isduplicate)
    fprintf(1,'All %u %s from %s already exist in %s---no new %s to add!\n',nits,thewhat,INPfile,dbname,thewhat);
    return
elseif sum(isduplicate) > 0
    if bevocal
        fprintf(1,'I found %u duplicate %s already in the database %s...?!\n',sum(isduplicate),thewhat,dbname)
        fprintf(1,'There are %u new %s to add to %s...\n',sum(~isduplicate),thewhat,dbname)
    end
end

% Select the maximum id already in the table
maxid = mysql_dbquery(dbc,sprintf('SELECT MAX(%s) FROM %s',theid,thetable));
maxid = maxid{1}; % the maximum id -- the new items will have ids greater than this
if isempty(maxid), maxid = 0; end

% Assemble and execute the INSERT queries
fprintf('Adding %u new %s to the %s table in %s...',sum(~isduplicate),thewhat,thetable,dbname)
switch importwhat
case 'ts' % Add time series to the TimeSeries table
    SQL_add_chunked(dbc,'INSERT INTO TimeSeries (FileName, Keywords, Length, Data) VALUES',toadd,isduplicate);
case 'ops' % Add operations to the Operations table
    SQL_add_chunked(dbc,'INSERT INTO Operations (OpName, Code, MasterLabel, Keywords) VALUES',toadd,isduplicate);        
case 'mops'
    SQL_add_chunked(dbc,'INSERT INTO MasterOperations (MasterLabel, MasterCode) VALUES',toadd,isduplicate);
end
fprintf(1,' Done!\n')

% Add new entries to the Results table where the TIMESERIES (ts_id) doesn't already exist
if ~strcmp(importwhat,'mops')
    resultstic = tic;
    if bevocal
        fprintf(1,'Updating the Results Table in %s (this could take a while, please be patient!)...',dbname)
    end
    switch importwhat
    case 'ts'
        [rs,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,m_id) SELECT t.ts_id,o.m_id FROM TimeSeries t' ...
                                ' CROSS JOIN Operations o ON t.ts_id > %u'],maxid));
    case 'ops'
        [rs,emsg] = mysql_dbexecute(dbc,sprintf(['INSERT INTO Results (ts_id,m_id) SELECT t.ts_id,o.m_id FROM TimeSeries t' ...
                                ' CROSS JOIN Operations o ON o.m_id > %u'],maxid));
    end
    if ~isempty(emsg),
        fprintf(1,' error. This is really not good.\n');
        fprintf(1,'%s\n',emsg);
        keyboard
    else
        if bevocal, fprintf(1,' initialized in %s!!\n',benrighttime(toc(resultstic))); end
    end
end

if ~strcmp(importwhat,'mops')
    % Update the keywords table
    fprintf(1,'Updating the %s table in %s...',thektable,dbname)

    % First find unique keywords from new time series by splitting against commas
    switch importwhat
    case 'ts'
        kws = {timeseries(~isduplicate).Keywords};
    case 'ops'
        kws = {operation(~isduplicate).Keywords};
    end
    
    kwsplit = cell(length(kws),1); % split into each individual keyword
    ukws = {};
    for i = 1:length(kws)
        kwsplit{i} = regexp(kws{i},',','split','ignorecase');
        for j = 1:length(kwsplit{i})
            if ~ismember(kwsplit{i}{j},ukws) % add it to ukws
                ukws{end+1} = kwsplit{i}{j};
            end
        end
    end
    nkw = length(ukws); % the number of unique keywords in the new set of time series
    if bevocal, fprintf(1,'\nI found %u unique keywords in the %s in %s...',nkw,thewhat,INPfile); end

    % How many overlap with existing keywords??
    allkws = mysql_dbquery(dbc,sprintf('SELECT Keyword FROM %s',thektable));
    if ~isempty(allkws) % the table may be empty, in which case all keywords will be new
        isnew = cellfun(@(x)~isempty(x),regexp(ukws,allkws,'ignorecase')); % ignore case for keywords
    else
        isnew = ones(nkw,1); % all are new
    end
    
    if sum(isnew) > 0
        if bevocal
            fprintf(1,['\nIt turns out that %u keywords are completely new and will be added ' ...
                        'to the %s table in %s'],sum(isnew),thektable,dbname)
        end
        % Add the new keywords to the Keywords table
        insertstring = sprintf('INSERT INTO %s (Keyword,NumOccur) VALUES',thektable);
        toadd = cell(sum(isnew),1);
        fisnew = find(isnew); % indicies of new keywords
        for k = 1:sum(isnew);
            toadd{k} = sprintf('(''%s'',0)',ukws{fisnew(k)});
        end
        SQL_add_chunked(dbc,insertstring,toadd);
        fprintf(1,' added %u new keywords!\n',nkw)
    else
        if bevocal
            fprintf(1,['\nIt turns out that all new keywords already exist in ' ...
                        'the %s table in %s -- there are no new keywords to add\n'],sum(isnew),thektable,dbname)
        end
    end
    
    % Fill new relationships
    fprintf(1,'Writing new keyword relationships to the %s table in %s...',thereltable,dbname)
    % Try doing it from scratch

    % allkws, allids
    switch importwhat
    case 'ts'
        allnames = bencat({timeseries.Filename},',','''');
    case 'ops'
        allnames = bencat({operation.Name},',','''');
    end
    ourids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE %s IN (%s)',theid,thetable,thename,allnames));
    ourids = vertcat(ourids{:}); % ids matching FileNames/OpNames
    ourkids = mysql_dbquery(dbc,sprintf('SELECT %s FROM %s WHERE Keyword IN (%s)',thekid,thektable,bencat(ukws,',','''')));
    ourkids = vertcat(ourkids{:}); % ids matching FileNames/OpNames
    nkwrels = sum(cellfun(@(x)length(x),kwsplit)); % number of keyword relationships in the input file
    addcell = cell(nkwrels,1);
    ii = 1;
    for i = 1:nits
        for j = 1:length(kwsplit{i})
            addcell{ii} = sprintf('(%u,%u)',ourids(i),ourkids(strcmp(kwsplit{i}{j},ukws)));
            ii = ii+1;
        end
    end
    SQL_add_chunked(dbc,sprintf('INSERT INTO %s (%s,%s) VALUES',thereltable,theid,thekid),addcell); % add them all in chunks
    
    % REALLY SLOW:
    % for i = 1:length(kws) % each time series or operation
    %     for j = 1:length(kwsplit{i}) % each keyword assinged to the time series or operation
    %         switch importwhat
    %         case 'ts'
    %             InsertString = sprintf(['INSERT INTO TsKeywordsRelate (tskw_id, ts_id) SELECT ' ...
    %                 'tskw_id, ts_id FROM TimeSeriesKeywords, TimeSeries ' ...
    %                 'WHERE TimeSeriesKeywords.Keyword = ''%s'' AND TimeSeries.Filename = ''%s'''], ...
    %                             kwsplit{i}{j},timeseries(i).Filename);
    %             [~,emsg] = mysql_dbexecute(dbc, InsertString);
    %             if ~isempty(emsg)
    %                 fprintf(1,'\nError inserting %s,%s to TsKeywordsRelate\n',...
    %                                 kwsplit{i}{j},timeseries(i).Filename); keyboard
    %             end
    %         case 'ops'
    %             InsertString = sprintf(['INSERT INTO OpKeywordsRelate (mkw_id, m_id) SELECT ' ...
    %                 'mkw_id, m_id FROM OperationKeywords, Operations ' ...
    %                 'WHERE OperationKeywords.Keyword = ''%s'' AND Operations.OpName = ''%s'''], ...
    %                             kwsplit{i}{j},esc(operation(i).Name));
    %             [~,emsg] = mysql_dbexecute(dbc, InsertString);
    %             if ~isempty(emsg)
    %                 fprintf(1,'\nError inserting %s,%s to OpKeywordsRelate\n',...
    %                             kwsplit{i}{j},operation(i).Name); keyboard
    %             end
    %         end
    %     end
    % end
    
    % Increment Nmatches in the keywords table
    fprintf(1,' done.\nNow the match counts...')
    % Redo them from scratch should be easier actually
    for k = 1:nkw % keywords implicated in this import
        SelectString = sprintf('(SELECT %s FROM %s WHERE Keyword = ''%s'')',thekid,thektable,ukws{k});
        themkw = mysql_dbquery(dbc,SelectString);
        UpdateString = sprintf('UPDATE %s SET NumOccur = (SELECT COUNT(*) FROM %s WHERE %s = %u) WHERE %s = %u', ...
                                    thektable,thereltable,thekid,themkw{1},thekid,themkw{1});
        [~,emsg] = mysql_dbexecute(dbc, UpdateString);
        if ~isempty(emsg)
            fprintf(1,'\n Error updating keyword count in %s',thektable)
            keyboard
        end
    end
    % for k = 1:nkw % for each unique keyword in the keyword table...
    %     % nnkw = sum(cellfun(@(x)ismember(ukws{k},x),kwsplit));
    %     Selectmkwid = sprintf('(SELECT %s FROM %s WHERE Keyword = ''%s'')',thekid,thektable,ukws{k});
    %     SelectCount = sprintf(['SELECT COUNT(*) FROM %s WHERE %s = %s ' ...
    %                             'AND %s > %u'],thereltable,thekid,Selectmkwid,theid,maxid);
    %     UpdateString = sprintf(['UPDATE %s SET NumOccur = NumOccur + (%s) ' ...
    %                             'WHERE Keyword = ''%s'''],thektable,SelectCount,ukws{k});
    %     [~,emsg] = mysql_dbexecute(dbc, UpdateString);
    %     if ~isempty(emsg)
    %         keyboard
    %         fprintf(1,'\n Error updating keyword count in %s',thektable)
    %     end
    % end
    fprintf(1,'done\n')
end

% Update master/operation links
if ismember(importwhat,{'mops','ops'}) % there may be new links
    % Delete the linking table and recreate it from scratch is easiest
    fprintf(1,'Filling MasterPointerRelate...');
    mysql_dbexecute(dbc,'DELETE FROM MasterPointerRelate');
    InsertString = ['INSERT INTO MasterPointerRelate select m.mop_id,o.m_id FROM MasterOperations ' ...
                            'm JOIN Operations o ON m.MasterLabel = o.MasterLabel'];
    [~,emsg] = mysql_dbexecute(dbc,InsertString);
    if isempty(emsg)
        fprintf(' done\n');
    else
        fprintf(1,' shit. Error joining the MasterOperations and Operations tables:\n');
        fprintf(1,'%s\n',emsg); keyboard
    end
    
    % if strcmp(importwhat,'ops')
    %     % operations were imported -- match their MasterLabels with elements of the MasterOperations table using mySQL JOIN
    %     InsertString = ['INSERT INTO MasterPointerRelate SELECT m.mop_id,o.m_id FROM MasterOperations m JOIN ' ...
    %                         'Operations o ON m.MasterLabel = o.MasterLabel WHERE o.m_id > %u',maxid];
    % else
    %     InsertString = ['INSERT INTO MasterPointerRelate SELECT m.mop_id,o.m_id FROM MasterOperations m JOIN ' ...
    %                         'Operations o ON m.MasterLabel = o.MasterLabel WHERE m.mop_id > %u',maxid];
    % end
    
    M_ids = mysql_dbquery(dbc,'SELECT mop_id FROM MasterOperations');
	M_ids = vertcat(M_ids{:}); % vector of master_ids    
    for k = 1:length(M_ids)
    	UpdateString = sprintf(['UPDATE MasterOperations SET NPointTo = ' ...
							'(SELECT COUNT(mop_id) FROM MasterPointerRelate WHERE mop_id = %u)' ...
        					'WHERE mop_id = %u'],M_ids(k),M_ids(k));
    	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);
    	if ~isempty(emsg)
    		fprintf(1,'Error counting NPointTo operations for mop_id = %u\n',M_ids(k)); keyboard
    	end
    end
end

%% Close database
SQL_closedatabase(dbc)

fprintf('All tasks completed reading %s for %s into %s in %s.\n',INPfile,thewhat,dbname,benrighttime(toc(ticker)));

end