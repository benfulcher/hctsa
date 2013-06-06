% function SQL_import_timeseries(whatfile, INPfile, dbname)
% INPUTS:
% whatfile -- the type of input file (time series, sources, categories, distributioncodes)
% INPfile -- the filename of the textfile to be read in [default = INP_ts.txt or INP_mets.txt]
% The input file should be formatted with spaces as delimiters between the entries, and in the order specified below
% the column labels of the TimeSeries or Operations table in the database are in this order:
% Ben Fulcher June 2013




esc = @sqlescapestring; % We will use this function alot to parse strings to and from the database

%% 1. Use mySQL commands to read in the data from local .csv files
% run import.sql

TablesToFill = {'TimeSeriesCategories', ... % TimeSeriesCategories Table
            'TimeSeriesDistributionCodes', ...    % Time Series distribution codes Table
            'TimeSeriesSource', ...     % TimeSeriesSource Table
            'TimeSeries'};           % TimeSeries Table

nTablesToFill = length(TablesToFill);
% Use SQL_ImportDataString to get mySQL code for importing all this data
ImportDataStrings = arrayfun(@(x)SQL_ImportDataString(TablesToFill{x}),1:nTablesToFill,'UniformOutput',0);

fprintf(1,'%d Tables to fill\n',nTablesToFill)
[dbc, dbname] = SQL_opendatabase; % opens dbc, the default database (named dbname)
fprintf(1,'Filling tables in %s from local .csv files\n',dbname);
for j = 1:nTablesToFill
    [rs,emsg] = mysql_dbexecute(dbc,ImportDataStrings{j});
    if ~isempty(rs)
        fprintf(1,'Imported data into table: %s\n',TablesToFill{j});
    else
        fprintf(1,'**** Error importing data into table: %s\n',TablesToFill{j});
        fprintf(1,'%s',emsg);
    end
end

fprintf(1,'Tables in %s filled successfully\n',dbname);


%% 2. Link Source_ids
% Go through each source and copy its id to the TimeSeries Table
SelectString = 'SELECT Source_id, SourceName FROM TimeSeriesSource';
[rs,~,~,emsg] = mysql_dbquery(dbc,SelectString);
Sourceids = vertcat(rs{:,1});
SourceNames = rs(:,2);
nSources = length(SourceNames);
fprintf(1,'We have %g sources\n',nSources);
for j = 1:nSources
    nentries = mysql_dbquery(dbc,['SELECT COUNT(*) FROM TimeSeries WHERE SourceString = ''' esc(SourceNames{j}) '''']);
    nentries = nentries{1};
    UpdateString = ['UPDATE TimeSeries SET Source_id = ' num2str(Sourceids(j)) ' WHERE SourceString = ''' esc(SourceNames{j}) ''''];
    [rs,emsg] = mysql_dbexecute(dbc,UpdateString);
    if isempty(emsg)
        fprintf(1,'Updated %g entries in TimeSeries by assigning %s to a Source_id %g \n',nentries,SourceNames{j},Sourceids(j))
    else
        fprintf(1,'Error updating entries in TimeSeries for %s\n',SourceNames{j})
    end
end
fprintf(1,'******All SourceStrings in TimeSeries Table Linked Successfully*********\n')

% Fill NMembers for Sources [[[NONONO INCORPORATE INTO ABOVE....]]]
for j = 1:nSources
    nentries = mysql_dbquery(dbc,['SELECT COUNT(*) FROM TimeSeries WHERE Source_id = ''' num2str(Sourceids(j)) '''']);
    nentries = nentries{1};
    UpdateString = ['UPDATE TimeSeriesSource SET NMembers = ' num2str(nentries) ' WHERE Source_id = ' num2str(Sourceids(j))];
    [rs,emsg] = mysql_dbexecute(dbc,UpdateString);
    if isempty(emsg)
        fprintf(1,'Updated Nmembers in TimeSeriesSource: %g for Source_id %g \n',nentries,Sourceids(j))
    else
        fprintf(1,'Error updating Nmembers in TimeSeriesSource for %s\n',SourceNames{j})
    end
end

%% 3. Link Category_ids
SelectString = 'SELECT Category_id, CategoryName FROM TimeSeriesCategories';
[rs,~,~,emsg] = mysql_dbquery(dbc,SelectString);
Categoryids = vertcat(rs{:,1});
CategoryNames = rs(:,2);
nCategories = length(CategoryNames);
fprintf(1,'We have %g categories\n',nCategories);
for j = 1:nSources
    [nentries,~,~,emsg] = mysql_dbquery(dbc,['SELECT COUNT(*) FROM TimeSeries WHERE CategoryString = ''' esc(CategoryNames{j}) '''']);
    nentries = nentries{1};
    UpdateString = ['UPDATE TimeSeries SET Category_id = ' num2str(Categoryids(j)) ' WHERE CategoryString = ''' esc(CategoryNames{j}) ''''];
    [rs,emsg] = mysql_dbexecute(dbc,UpdateString);
    if isempty(emsg)
        fprintf(1,'Updated %g entries in TimeSeries by assigning %s to a Category_id %g \n',nentries,CategoryNames{j},Categoryids(j))
    else
        fprintf(1,'Error updating entries in TimeSeries for %s\n',CategoryNames{j})
    end
end
fprintf(1,'****All CategoryStrings in TimeSeries Table Linked Successfully\n')


SQL_closedatabase(dbc) % close the connection to the database


% tic;
% 
% % whatfile
% if nargin < 1
%     error('No input specified: what type of file am I looking for??')
% end
% validinps = {'ts','source','category','codes'};
% if ~ismember(whatfile,validinps)
%     error('Invalid input specifier')
% end
% 
% % INPfile
% if nargin < 2 || isempty(INPfile)
%     INPfile = 'INP_ts.txt';
% end
% disp(['Using input file ' INPfile]);
% 
% % dbname
% if nargin < 3
%     dbname = []; % uses default database, as specified in SQL_opendatabase
% end
% 
% %% Open Database
% [dbc, dbname] = SQL_opendatabase(dbname);


% %% Open and read file
% fid = fopen(INPfile); % opens the file using fopen
% 
% switch whatfile
% case 'ts'
%     disp(['Reading ' INPfile ' for time series']);
%     disp('Format: Filename, Keywords, Quantity, Unit, Sampling Rate, Description, Source Name, CategoryName');
% 
%     % Read in data from file:
%     datain = textscan(fid,'%s %s %s %s %s %s %s %s','Delimiter','\t');
% 
%     % Store results in Nx1 cells of strings, where N is the number of inputs to be added
%     % Ignore the first row?
%     % daatin{1}; % filename strings
%     % datain{2}; % keywords
%     % datain{3}; % quantity measured
%     % datain{4}; % unit of the quantity measured
%     % datain{5}; % sampling rate of measurement
%     % datain{6}; % description of this particular time series
%     % datain{7}; % name of the (synthetic or real-world) source of this time series that should match an existing source name
%     % datain{8}; % Category to assign to the time series that should exist in the categories table
%     
%     disp('What do you think of this???:');
%     nlinesshow = 5; % show first and last nlinesshow
%     for i = 1:nlines
%         fprintf(1,'%s, %s, %s, %s, %s, %s, %s, %s\n',x{1}{i},x{2}{i},x{3}{i},x{4}{i},x{5}{i},x{6}{i},x{7}{i},x{8}{i})
%     end
%     fprintf(1,'............................\n')
%     for i = nlines:-1:0
%         fprintf(1,'%s, %s, %s, %s, %s, %s, %s, %s\n',x{1}{end-i},x{2}{end-i},x{3}{end-i},x{4}{end-i},x{5}{end-i},x{6}{end-i},x{7}{end-i},x{8}{end-i})
%     end
%     reply = input('Is the first row a header?? -- ''y'' to delete...','s');
%     if strcmp(reply,'y')
%         for i = 1:length(datain), datain{i}(1) = []; end % remove the first entry in each
%         for i = 1:nlinesshow
%             fprintf(1,'%s, %s, %s, %s, %s, %s, %s, %s\n',datain{1}{i},datain{2}{i},datain{3}{i},datain{4}{i},datain{5}{i},datain{6}{i},datain{7}{i},datain{8}{i})
%         end
%         fprintf(1,'............................\n')
%         for i = nlinesshow:-1:0
%             fprintf(1,'%s, %s, %s, %s, %s, %s, %s, %s\n',datain{1}{end-i},datain{2}{end-i},datain{3}{end-i},datain{4}{end-i},datain{5}{end-i},datain{6}{end-i},datain{7}{end-i},datain{8}{end-i})
%         end
%         input('Deleted, good enough!!?');
%     end
% 
%     
% case 'source'
%     disp(['Reading ' INPfile ' for time series sources']);
%     disp('Format: Source Name, Source Link, Source Description, Processing Notes, Added By, Distribution Code id, Parent Name');
% 
%     % Read in data from file:
%     datain = textscan(fid,'%s %s %s %s %s %s %s','Delimiter',',');
%     
% case 'codes'
%     disp(['Reading ' INPfile ' for time series sources']);
%     disp('Format: Code, Description');
% 
%     % Read in data from file:
%     datain = textscan(fid,'%f %s','Delimiter',',')
%     
%     for i = 1:nlinesshow
%         fprintf(1,'%f, %s\n',datain{1}{i},datain{2}{i},datain{3}{i},datain{4}{i},datain{5}{i},datain{6}{i},datain{7}{i},datain{8}{i})
%     end
% end
% 
% fclose(fid); % closes the file

% %% 
% 
% chunksize = 1000; % Maximum of 1000 entries added per DB operation
% esc = @sqlescapestring;
% 
% switch metorts
%     case 'ts'
%         %error('For now, use TSQ_add() for timeseries');
%         if length(tsf) == 0
%             error('Index file seems to be empty...');
%         end
%         toadd = cell(length(tsf),1);
%         fprintf('Indexing TimeSeries');
%         for j=1:length(tsf)
%             timeseries(j).Filename = tsf{j};
%             timeseries(j).Keywords = tskw{j};
%             x = dlmread(tsf{j});
%             timeseries(j).Length = length(x);
%             toadd{j} = sprintf('(''%s'', ''%s'',%i, NOW())',esc(timeseries(j).Filename),esc(timeseries(j).Keywords),timeseries(j).Length);
%         end
%         
%         existing = mysql_dbquery(dbc,'Select Filename from TimeSeries');
%         isduplicate = ismember({timeseries.Filename},existing); % isduplicate is 1 if the MasterLabel already exists
% 
%         % Assemble the query
%         if all(isduplicate)
%             fprintf(', no new TimeSeries to add\n');
%         else
%             fprintf(', adding %i new entries...',sum(~isduplicate));
%             SQL_add_chunked(dbc,'INSERT INTO TimeSeries (FileName, Keywords, Length, LastModified) VALUES',toadd,isduplicate);
%             fprintf('done!\n');
%             
%             % Add new entries to the Results table where the TIMESERIES (ts_id) doesn't already exist  
%             fprintf('Updating the Results table (this could take a while, BE PATIENT) ...')
%             for j= 1:length(tsf)
%                 if ~isduplicate(j)
%                     mysql_dbexecute(dbc,sprintf('INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o WHERE t.FileName=''%s''',timeseries(j).Filename));
%                 end
%             end
%             % This query needs to be redone when the table is huge
%             %mysql_dbexecute(dbc,'INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o LEFT JOIN Results r ON r.ts_id=t.ts_id WHERE r.ts_id IS NULL ORDER BY ts_id,m_id');
%             fprintf('done!\n')
%             
%             fprintf('Updating the Timeseries Keywords table...')
%             SQL_update_tskw(dbname)
%             fprintf('done!\n')
%         end
%         
%         if any(isduplicate)
%             warning('Duplicate timeseries were found. If the keywords have changed, you will need to run TSQ_add() to update them');
%         end
%         
%     case 'mets'
%         
%         update_secondary = 0;
%         
%         %%% METRICS
%         %% Still more formatting to do: Masters first
%         % Masters
%         Mr = strmatch('M', methin); % those lines that start with this string pattern
%         nmM = length(Mr); % number of master rows
%         
%         if nmM > 0
%             fprintf('Indexing MasterOperations');
%             master_text = methin(Mr);
%             methin(Mr) = ''; % remove master rows from methin
%             toadd = cell(size(master_text,1),1);
%             for j = 1:size(master_text,1)
%                 tmp = textscan(master_text{j},'%s %s %s',1);
%                 if any(cellfun(@isempty, tmp))
%                     fprintf('%s is incorrectly formatted (read <%s>, <%s>, <%s>',master_text{j},tmp{1},tmp{2},tmp{3});
%                     beep;
%                     keyboard;
%                 end
%                 master(j).type = tmp{1}{1};
%                 master(j).MasterCode = tmp{2}{1};
%                 master(j).MasterLabel = tmp{3}{1};
%                 toadd{j} = sprintf('(''%s'', ''%s'', NOW())',esc(master(j).MasterLabel),esc(master(j).MasterCode));
%             end
% 
%             % Now we need to check for duplicates
%             existing = mysql_dbquery(dbc,'SELECT MasterLabel FROM MasterOperations');
%             isduplicate = ismember({master.MasterLabel},existing); % isduplicate is 1 if the MasterLabel already exists
%             
%             % Assemble the query
%             if all(isduplicate)
%                 fprintf(', no new MasterOperations to add\n');
%             else
%                 fprintf(', adding %i new entries...',sum(~isduplicate));
%                 SQL_add_chunked(dbc,'INSERT INTO MasterOperations (MasterLabel, MasterCode, LastModified) VALUES',toadd,isduplicate);
%                 fprintf('done!\n');
%                 update_secondary = 1;
%             end
%         else
%             fprintf('No MasterOperations were found\n');
%         end
% 
%         % Operations
%         ops = [strmatch('P', methin); strmatch('S', methin)]; % those lines that start with this string pattern
%         nops = length(ops); % number of master rows
% 
%         if nops>0
%             fprintf('Indexing Operations');            
%             op_text = methin(ops);
%             toadd = cell(size(op_text,1),1);
%             op_names = toadd;
%             
%             for j = 1:size(op_text,1)
%                 tmp = textscan(op_text{j},'%s %s %s %s %s',1);
%                 if any(cellfun(@isempty, tmp))
%                     fprintf('%s is incorrectly formatted',op_text{j},tmp{1},tmp{2},tmp{3},tmp{4});
%                     beep;
%                     keyboard;
%                 end
%                 Pointer = strcmp(tmp{1}{1},'P');
%                 OpCode = tmp{2}{1};
%                 op_names{j} = tmp{3}{1};
%                 %op(j).MasterLabel = strtok(op(j).OpCode,'.');
%                 MasterLabel = OpCode(1:strfind(OpCode,'.')-1);
%                 Normalize = 1;
%                 Keywords = tmp{5}{1}; % keywords 
%                 toadd{j} = sprintf('(''%s'', ''%s'',''%s'',%d,%d,''%s'', NOW())',esc(op_names{j}),esc(OpCode),esc(MasterLabel),Normalize,Pointer,esc(Keywords));
%             end
% 
%             % Now we need to check for duplicates
%             existing = mysql_dbquery(dbc,'SELECT OpName FROM Operations');
%             isduplicate = ismember(op_names,existing); % isduplicate is 1 if the MasterLabel already exists
% 
%             
%             % Assemble the query
%             if all(isduplicate)
%                 fprintf(', no new Operations to add\n');
%             else
%                 fprintf(', adding %i new entries...',sum(~isduplicate));
%                 SQL_add_chunked(dbc,'INSERT INTO Operations (OpName, Code, MasterLabel, Normalize, Pointer, Keywords, LastModified) VALUES',toadd,isduplicate);
%                 fprintf('done!\n');
%                 update_secondary = 1;
%             end
%         else
%             fprintf('No Operations were found\n');
%         end
% 
%         if update_secondary
%             fprintf('Updating the list of keywords...');
%             SQL_update_mkw(dbname) % update operation keywords
%             fprintf('done!\n');
%             
%             % Repopulate the MasterPointerRelate table
%             fprintf('Updating MasterPointerRelate...');
%             mysql_dbexecute(dbc,'DELETE FROM MasterPointerRelate;');
%             mysql_dbexecute(dbc,'INSERT INTO MasterPointerRelate select m.Master_id,o.m_id FROM MasterOperations m JOIN Operations o ON m.MasterLabel= o.MasterLabel;');
%             SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
%             fprintf('done!\n');
%             
%             % Add new entries to the Results table where the OPERATION (m_id) doesn't already exist  
%             fprintf('Updating the Results table (this could take a while, BE PATIENT) ...')
%             mysql_dbexecute(dbc,'INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o LEFT JOIN Results r ON r.m_id=o.m_id WHERE r.m_id IS NULL ORDER BY ts_id,m_id');    
%             fprintf('done!\n');
%         end
% end
% 
% %% Close database
% SQL_closedatabase(dbc)
% 
% fprintf('All tasks completed in %.0f seconds\n',toc);



% end