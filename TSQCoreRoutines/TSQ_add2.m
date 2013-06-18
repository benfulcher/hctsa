function TSQ_add2(metorts, INPfile, dbname)
% Compared to TSQ_add2
% This is UNABLE to UPDATE KEYWORDS
%
% So use normal TSQ_add for that

%%% TSQ_add
% Adds a time series or metric to the database with the specified columns, collabs, of the table given
% [MAYBE] NOT intended to cover changing entries in the database -- TSQ_alter should be used for this.
% 
% INPUTS: metorts -- 'mets' or 'ts' -- write to which of these?
% 		  inpfilename -- the filename of the textfile to be read in [default = INP_ts.txt or INP_mets.txt]
% 
% The input file should be formatted with spaces as delimiters between the entries, and in the order specified below
% the column labels of the TimeSeries or Operations table in the database are in this order:
% 			 'FileName','Keywords' for TimeSeries
% 			 'Pointer'(M,S,P),'Code','OpName','Keywords' for Operations (without keywords for master functions)
% 
% Ben Fulcher 3/12/2009
% Ben Fulcher 12/1/2010: added dbname option

%% Check Inputs, Set defaults
%% ROMESH- For now, add all of the mets first to avoid poor performance
tic;

% metorts
if nargin < 1 || isempty(metorts) || ~(strcmp(metorts,'mets') || strcmp(metorts, 'ts'))
	disp('Error setting first input argument -- should be ''ts'' or ''mets''');
	return
end

% inpfilename
if nargin < 2 || isempty(INPfile)
	if strcmp(metorts,'ts')
		INPfile = 'INP_ts.txt';
	else
		INPfile = 'INP_mets.txt';
	end
end
disp(['Using input file ' INPfile]);

% dbname
if nargin < 3
	dbname = []; % uses default database, as specified in SQL_opendatabase
end

%% Open Database
[dbc dbname] = SQL_opendatabase(dbname);


%% Open, read file
fid = fopen(INPfile);

if strcmp(metorts,'ts')
	datin = textscan(fid,'%s %s','CommentStyle','%');
	tsf = datin{1}; % filename strings -- Nx1 cell of strings
	tskw = datin{2}; % keywords -- Nx1 cell of string cells
	if any(cellfun(@isempty,tsf)) || any(cellfun(@isempty,tskw))
		disp('NO NO NO!!: badly formatted...'); beep
		keyboard;
    end
else
	methin = textscan(fid,'%s','Delimiter','','CommentStyle','%');
	methin = methin{1};	
end
fclose(fid);

%% 

chunksize = 1000; % Maximum of 1000 entries added per DB operation
esc = @sqlescapestring;

switch metorts
    case 'ts'
        %error('For now, use TSQ_add() for timeseries');
        if length(tsf) == 0
            error('Index file seems to be empty...');
        end
        toadd = cell(length(tsf),1);
        fprintf('Indexing TimeSeries');
        for j=1:length(tsf)
            timeseries(j).Filename = tsf{j};
            timeseries(j).Keywords = tskw{j};
            x = dlmread(tsf{j});
            timeseries(j).Length = length(x);
            toadd{j} = sprintf('(''%s'', ''%s'',%i, NOW())',esc(timeseries(j).Filename),esc(timeseries(j).Keywords),timeseries(j).Length);
        end
        
        existing = mysql_dbquery(dbc,'Select Filename from TimeSeries');
        isduplicate = ismember({timeseries.Filename},existing); % isduplicate is 1 if the MasterLabel already exists

        % Assemble the query
        if all(isduplicate)
            fprintf(', no new TimeSeries to add\n');
        else
            fprintf(', adding %i new entries...',sum(~isduplicate));
            SQL_add_chunked(dbc,'INSERT INTO TimeSeries (FileName, Keywords, Length, LastModified) VALUES',toadd,isduplicate);
            fprintf('done!\n');
            
            % Add new entries to the Results table where the TIMESERIES (ts_id) doesn't already exist  
            fprintf('Updating the Results table (this could take a while, BE PATIENT) ...')
            for j= 1:length(tsf)
                if ~isduplicate(j)
                    mysql_dbexecute(dbc,sprintf('INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o WHERE t.FileName=''%s''',timeseries(j).Filename));
                end
            end
            % This query needs to be redone when the table is huge
            %mysql_dbexecute(dbc,'INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o LEFT JOIN Results r ON r.ts_id=t.ts_id WHERE r.ts_id IS NULL ORDER BY ts_id,m_id');
            fprintf('done!\n')
            
            fprintf('Updating the Timeseries Keywords table...')
            SQL_update_tskw(dbname)
            fprintf('done!\n')
        end
        
        if any(isduplicate)
            warning('Duplicate timeseries were found. If the keywords have changed, you will need to run TSQ_add() to update them');
        end
        
    case 'mets'
        
        update_secondary = 1;
        
        %%% METRICS
        %% Still more formatting to do: Masters first
        % Masters
        Mr = strmatch('M', methin); % those lines that start with this string pattern
        nmM = length(Mr); % number of master rows
        
        if nmM>0
            fprintf('Indexing MasterOperations');
            master_text = methin(Mr);
            methin(Mr) = ''; % remove master rows from methin
            toadd = cell(size(master_text,1),1);
            for j = 1:size(master_text,1)
                tmp = textscan(master_text{j},'%s %s %s',1);
                if any(cellfun(@isempty, tmp))
                    fprintf('%s is incorrectly formatted (read <%s>, <%s>, <%s>',master_text{j},tmp{1},tmp{2},tmp{3});
                    beep;
                    keyboard;
                end
                master(j).type = tmp{1}{1};
                master(j).MasterCode = tmp{2}{1};
                master(j).MasterLabel = tmp{3}{1};
                toadd{j} = sprintf('(''%s'', ''%s'', NOW())',esc(master(j).MasterLabel),esc(master(j).MasterCode));
            end

            % Now we need to check for duplicates
            existing = mysql_dbquery(dbc,'SELECT MasterLabel FROM MasterOperations');
            isduplicate = ismember({master.MasterLabel},existing); % isduplicate is 1 if the MasterLabel already exists
            
            % Assemble the query
            if all(isduplicate)
                fprintf(', no new MasterOperations to add\n');
            else
                fprintf(', adding %i new entries...',sum(~isduplicate));
                SQL_add_chunked(dbc,'INSERT INTO MasterOperations (MasterLabel, MasterCode, LastModified) VALUES',toadd,isduplicate);
                fprintf('done!\n');
                update_secondary = 1;
            end
        else
            fprintf('No MasterOperations were found\n');
        end

        % Operations
        ops = [strmatch('P', methin); strmatch('S', methin)]; % those lines that start with this string pattern
        nops = length(ops); % number of master rows

        if nops>0
            fprintf('Indexing Operations');            
            op_text = methin(ops);
            toadd = cell(size(op_text,1),1);
            op_names = toadd;
            
            for j = 1:size(op_text,1)
                tmp = textscan(op_text{j},'%s %s %s %s %s',1);
                if any(cellfun(@isempty, tmp))
                    fprintf('%s is incorrectly formatted',op_text{j},tmp{1},tmp{2},tmp{3},tmp{4});
                    beep;
                    keyboard;
                end
                Pointer = strcmp(tmp{1}{1},'P');
                OpCode = tmp{2}{1};
                op_names{j} = tmp{3}{1};
                %op(j).MasterLabel = strtok(op(j).OpCode,'.');
                MasterLabel = OpCode(1:strfind(OpCode,'.')-1);
                Keywords = tmp{5}{1}; % keywords 
                toadd{j} = sprintf('(''%s'', ''%s'',''%s'',%d,''%s'')',esc(op_names{j}),esc(OpCode),esc(MasterLabel),Pointer,esc(Keywords));
            end

            % Now we need to check for duplicates
            existing = mysql_dbquery(dbc,'SELECT OpName FROM Operations');
            isduplicate = ismember(op_names,existing); % isduplicate is 1 if the MasterLabel already exists

            
            % Assemble the query
            if all(isduplicate)
                fprintf(', no new Operations to add\n');
            else
                fprintf(', adding %i new entries...',sum(~isduplicate));
                SQL_add_chunked(dbc,'INSERT INTO Operations (OpName, Code, MasterLabel, Pointer, Keywords) VALUES',toadd,isduplicate);
                fprintf('done!\n');
                update_secondary = 1;
            end
        else
            fprintf('No Operations were found\n');
        end

        if update_secondary
            fprintf('Updating the list of keywords...');
            SQL_update_mkw(dbname) % update operation keywords
            fprintf('done!\n');
            
            % Repopulate the MasterPointerRelate table
            fprintf('Updating MasterPointerRelate...');
            mysql_dbexecute(dbc,'DELETE FROM MasterPointerRelate;');
            mysql_dbexecute(dbc,'INSERT INTO MasterPointerRelate select m.mop_id,o.m_id FROM MasterOperations m JOIN Operations o ON m.MasterLabel = o.MasterLabel;');
            SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
            fprintf('done!\n');
            
            % Add new entries to the Results table where the OPERATION (m_id) doesn't already exist  
            fprintf('Updating the Results table (this could take a while, BE PATIENT) ...')
            mysql_dbexecute(dbc,'INSERT INTO Results (ts_id,m_id)  SELECT t.ts_id,o.m_id FROM TimeSeries t CROSS JOIN Operations o LEFT JOIN Results r ON r.m_id=o.m_id WHERE r.m_id IS NULL ORDER BY ts_id,m_id');	
            fprintf('done!\n');
        end
end

%% Close database
SQL_closedatabase(dbc)

fprintf('All tasks completed in %.0f seconds\n',toc);
