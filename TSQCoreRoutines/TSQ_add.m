function TSQ_add(metorts, INPfile, dbname)

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

% metorts
if nargin < 1 || isempty(metorts) || ~(strcmp(metorts,'mets') || strcmp(metorts, 'ts'))
	disp(['Error setting first input argument -- should be ''ts'' or ''mets''']);
	return
end
disp(['Using ' metorts]);

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
if nargin<3
	dbname = []; % uses default database, as specified in SQL_opendatabase
end

% Decided against this input -- too difficult
% % collabs, column labels to write to
% if nargin<2 || isempty(collabs)
% 	if strcmp(metorts,'ts')
% 		collabs = {'FileName','Keywords'}
% 	else
% 		collabs = {'OpName','Pointer','Keywords','Code'};
% 	end
% end
% disp(['Inserting into columns ' bencat(collabs)]);


%% Open Database
[dbc dbname] = SQL_opendatabase(dbname);

input([dbname ' ok? (CONTROL-C NOW IF NOT!!)']);

%% Open, read file
fid = fopen(INPfile);

if strcmp(metorts,'ts')
	datin = textscan(fid,'%s %s','CommentStyle','%');
	% datin = textscan(fid,'%d %s %s','CommentStyle','%');
	% tsid = datin{1};
	tsf = datin{1}; % filename strings -- Nx1 cell of strings
	tskw = datin{2}; % keywords -- Nx1 cell of string cells
	% tspp = datin{3}; % preprocessings
	if any(cellfun(@isempty,tsf)) || any(cellfun(@isempty,tskw))
		disp('NO NO NO!!: badly formatted...'); beep
		keyboard;
    end
else
	methin = textscan(fid,'%s','Delimiter','','CommentStyle','%');
	methin = methin{1};	
end
fclose(fid);


disp(['Successfully read from the input file: ' INPfile]);

%% 
addcount = 0;

if strcmp(metorts, 'ts')
	
	%% Read input text file and map to database

	nts_add = length(tsf);
	ts_ids = zeros(nts_add,1);
    
	for i=1:nts_add
		SelectString = ['SELECT Keywords FROM TimeSeries WHERE Filename = ''' tsf{i} ''''];
		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc, SelectString);

		if isempty(qrc) && isempty(emsg) % no matches -- this is in INP_ts.txt but not in the mySQL database
			% disp('ABOUT TO INSERT A NEW ONE -- I KNOW FOR A FACT YOU DON''T WANT THIS!!')
			% keyboard
			
			%% Insert the new entry into the TimeSeries Table
			insertstring = ['INSERT INTO TimeSeries (FileName, Keywords, LastModified) VALUES ('''...
			 				tsf{i} ''', ''' tskw{i} ''', NOW() )'];
			% insertstring = ['INSERT INTO TimeSeries (FileName, Keywords, Preprocess, LastModified) VALUES ('''...
			%  				tsf{i} ''', ''' tskw{i} ''', ' num2str(tspp{i}) ', NOW() )'];
		    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
			if ~isempty(emsg)
				disp(['Error adding time series ' tsf{i} ' to database']); keyboard
			else
				%% Add new entries in Results Table
				% Get the ts_id of the new TimeSeries
				SelectString = ['SELECT ts_id from TimeSeries WHERE Filename = ''' tsf{i} ''''];
				[theid,qrf,rs,emsg] = mysql_dbquery(dbc, SelectString);
				theid = theid{1};
                ts_ids(i) = theid;
				
				insertstring = ['INSERT INTO Results (ts_id, m_id) SELECT ' num2str(theid) ', m_id FROM Operations'];
			    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
				if isempty(emsg)
					addcount = addcount + 1; % successfully added this time series
					disp(['Added ' tsf{i} ' (' num2str(i) '/' num2str(nts_add) ') to TimeSeries and initialized Results']);
				else
					disp(['Error initializing results for ' tsf{i} ' -- SHIT!!!: We''ve added the TimeSeries but the Results table is now inconsistent...']);
					disp(['Should really delete ' tsf{i}]);
					keyboard
				end
			end
		else
			if isempty(emsg) % no match, no error -- already exists in database
				input([tsf{i} ' already exists in the mySQL database. Checking up to date...']);
				
				if ~strcmp(tskw{i},qrc{1})
					% time series keywords changed from qrc{1} (in database) to tskw{1} (in input file)
					fprintf('%s \n\t\t %s \n %s \n\t\t %s \n','***** Keywords changed -- should I overwrite:',qrc{1},'with:',tskw{i});
					reply = input([' for ' tsf{i} ' ??'],'s');
					if ~strcmp(reply,'n')
						updatestring = ['UPDATE TimeSeries SET Keywords = ''' tskw{i} ''' WHERE FileName = ''' tsf{i} ''''];
					    [rs,emsg] = mysql_dbexecute(dbc, updatestring);
						if ~isempty(emsg)
							disp('error updating...'); keyboard
						end
					end
				end
				
			else % error
				input(['Error retrieving ' tsf{i}]); keyboard
			end
		end
	end


	if addcount>0
		disp(['Added ' num2str(addcount) ' ( / ' num2str(nts_add) ' ) new time series to database successfully! Really really!']);
		
		% Update:
		% (1) Metadata (length, positive-only,)
		disp('Writing Metadata');
		SQL_writetsmeta(dbname)
		
		% (2) Keywords
		disp('Updating TimeSeries keywords'); tic
		SQL_update_tskw(dbname)
		disp(['TimeSeries keywords updated (' benrighttime(toc) ')']);
		
		% (3) Percentage Calculated
		reply = input('Shall I update the percentage calculated and related columns? This takes forever, don''t do it! (''y'' for yes)','s');
		if strcmp(reply,'y')
			tic
			SQL_fillfromresults(ts_ids,[],[1 1 1],[1 1 0 0],dbname)
			disp(['Percentage calculated results took ' benrighttime(toc) ' to update']);
		end
		
		disp('Done and dusted!! Yeah!');
		
	else
		disp(['No time series out of the attempted ' num2str(nts_add) ' were actually added to the database :-( Something''s wrong....']);
	end
	
	
	
else
	%%% METRICS
	%% Still more formatting to do: Masters first
	% Masters
	Mr = strmatch('M', methin); % those lines that start with this string pattern
	nmM = length(Mr); % number of master rows
	
	if nmM>0
		masters = methin(Mr);
		char(masters)
		input(['There are ' num2str(nmM) ' master functions in the input file. Adding them now with your permission...']);
		methin(Mr) = ''; % remove master rows from methin
		masteraddcount = 0;
		% write them
		disp(['Formatting Master Functions First!!!']);
		for i=1:nmM
			mastertmp = textscan(masters{i},'%s %s %s',1);
			if any(cellfun(@isempty, mastertmp))
				disp(['NO NO NO!!: ' num2str(i) ' is a badly formatted master metric...']); beep
				keyboard;
		    end
			mastercode = sqlescapestring(char(mastertmp{2}));
			masterlabel = sqlescapestring(char(mastertmp{3}));

			matchcheckstring = ['SELECT MasterLabel, MasterCode FROM MasterOperations WHERE MasterLabel = ''' masterlabel ''''];
			[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,matchcheckstring);

			if isempty(qrc) % a new one -- insert it
				insertstring = ['INSERT INTO MasterOperations ' ...
								'(MasterLabel, MasterCode, LastModified) VALUES ' ...
								'(''' masterlabel ''', ''' mastercode ''', NOW() )'];
			    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
				if isempty(rs) % failed
					disp(['Inserting' masterlabel ' failed :(']); keyboard
				else
				    % successfully inserted new master metric
					masteraddcount = masteraddcount + 1;
				end
			else
				disp(['Already an entry for ' masterlabel]);
				if ~strmatch(qrc{2},mastercode) % this masterlabel already exists: but code doesn't match
					input(['Existing entry for ' masterlabel ' has inconsistent code! Ok?']);
				end
			end
		end
	
		disp(['Just added ' num2str(masteraddcount) ' ( / ' num2str(nmM) ') master functions to the database']);
	else
		disp(['No Master functions in the input file. Getting on with it.']);
	end
	
	nm = length(methin); % number of metrics (scalar integer)
	opaddcount = 0; % number of operations successfully added
	opkwchangecount = 0; % number of operations whose keywords are changed
	input(['There are ' num2str(nm) ' metrics in the input file... Let''s go one-by-one, shall we?'])
	tic
	for i=1:nm
		%% Format/retrieve relevant entry in INP_mets.txt
		methtmp = textscan(methin{i},'%s %s %s %n %s',1); % Delimiter is white space
		if any(cellfun(@isempty,methtmp))
			disp(['NO NO NO!!: metric ' num2str(i) ' poorly formatted...']); beep
			keyboard;
	    end

		% Pointer
		if strcmp(methtmp{1},'P')
	        oppoint=1; % a pointer
	    elseif strcmp(methtmp{1},'S')
			oppoint=0; % a single metric
		else
			disp(['No No No -- Error with pointer/single metric formatting']);
		end

		opcode = sqlescapestring(char(methtmp{2})); % code to call
		opname = char(methtmp{3}); % labels to associate
		opnormalize = methtmp{4}; % how to post-process/normalize
		opkws = char(methtmp{5}); % keywords

		%% Write to database
		matchcheckstring = ['SELECT OpName, Code, Normalize, Pointer, Keywords FROM Operations WHERE OpName = ''' opname ''''];
		[qrc,qrf,rs,emsg] = mysql_dbquery(dbc, matchcheckstring);

		if isempty(qrc) % a new one -- INSERT it
			insertstring = ['INSERT INTO Operations ' ...
							'(OpName, Code, Normalize, Pointer, Keywords, LastModified) VALUES (''' ...
			 				opname ''', ''' opcode ''', ' num2str(opnormalize) ', ' num2str(oppoint) ', ''' opkws ''', NOW() )'];
		    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
			if isempty(emsg)
				%% Add new entries to Results Table
				% Get the m_id of the new TimeSeries
				SelectString = ['SELECT m_id from Operations WHERE OpName = ''' opname ''''];
				[theid,qrf,rs,emsg] = mysql_dbquery(dbc, SelectString);
				theid = theid{1};
				insertstring = ['INSERT INTO Results (ts_id, m_id) SELECT ts_id, ' num2str(theid) ' FROM TimeSeries'];
			    [rs,emsg] = mysql_dbexecute(dbc, insertstring);
			
				if isempty(emsg)
					opaddcount = opaddcount + 1;
					disp(['Added ' opname ' to Operations and initialized Results (' num2str(i) '/' num2str(nm) ')']);
				else
					disp(['Error initializing results for ' opname ' -- SHIT!!!: We''ve added the Operation but the Results table is now inconsistent......']);
					disp(['Should really delete ' opname]);
					keyboard
				end				
			else
				disp(['Error inserting ' opname]); keyboard
			end
		else
			disp([opname ' already exists in the database... Checking up to date...']);
			% entry for this operation already exists -- check that information is still up to date

			% (1) Code changed -- kick up a fuss -- will affect existing results in storage
			if ~strmatch(qrc{2},opcode)
				disp(['CODE CHANGED FOR ' opname '!!!! I''m not touching this one: should be a new metric.']);
				disp('Skipping'); continue
			end

			% (2) Normalization changed -- no big deal, just update
			if qrc{3} ~= opnormalize
				% keyboard
				reply = input(['Normalization for ' opname ' changed from ' num2str(qrc{3}) ' to ' num2str(opnormalize) '. ' ...
						'Shall I update the database...? [n for ''no'']'],'s');
				if ~strcmp(reply,'n')
					updatestring = ['UPDATE Operations SET Normalize = ' num2str(opnormalize) ...
									', LastModified = NOW() WHERE OpName = ''' opname ''''];
				    [rs,emsg] = mysql_dbexecute(dbc, updatestring);
					if isempty(rs)
						disp(['Error updating normalization for ' opname]); keyboard
					else
						disp(['Updated ' opname ', which had an outdated normalization setting stored']);
					end
				else
					disp('Skipped');
				end
			end

			% (3) Pointer changed -- this should never happen
			if qrc{4}~=oppoint
				disp(['Pointer identity inconsistent with that in the database for ' opname '. This is weird...']); keyboard
			end

			% (4) Keywords changed -- ask to update
			if ~strcmp(qrc{5},opkws)
				reply = input(['Keywords changed for existing entry ' opname '. Shall I update the database...? [n for ''no'']'],'s');
				if ~strcmp(reply,'n')
					updatestring = ['UPDATE Operations SET Keywords = ''' opkws ...
									''', LastModified = NOW() WHERE OpName = ''' opname ''''];
				    [rs,emsg] = mysql_dbexecute(dbc, updatestring);
					if isempty(rs)
						disp(['Error updating keywords']); keyboard
					else
						disp(['Upated keywords for ' opname]);
						opkwchangecount = opkwchangecount + 1;
					end
				else
					disp(['Skipped']);
				end
			end
			
			
			
		end
	end
	
	if opaddcount>0
		disp(['Successfully added ' num2str(opaddcount) ' operations to the Operations table of the databae']);
		disp(['It took ' benrighttime(toc)]);

		% Update operation keywords
		disp('Doing Operation housekeeping'); tic
		SQL_update_mkw(dbname) % update operation keywords
		if nmM>0
			disp(['There were masters added -- updating linkage information']);
			SQL_linkpointermaster(dbname) % update master/pointer links
			SQL_masternpointto(dbname) % counts master/pointer links for MasterOperations table
		end
		disp(['Operation housekeeping took ' benrighttime(toc)]);

		
		% Reorder for optimization?
		% [should only added operations, since time series should have been added in order and with ordered m_ids]
		reply = input('The Results table may not be ordered... Want to reorder (IS EXTREMELY time consuming, maybe up to a week!... ''y'' to do this...)','s');
		if strcmp(reply,'y')
			disp('We''re doing it!!!');
			tic
			alterstring = 'ALTER TABLE Results ORDER BY ts_id, m_id';
		    [rs,emsg] = mysql_dbexecute(dbc, alterstring);
			if ~isempty(emsg)
				disp('Error reordering Results Table');
			else
				disp(['We did it!!! It only took ' benrighttime(toc)]);
			end
		end

	else
		disp(['None of the ' num2str(nm) ' operations were actually added :-(']);
		
	end
	
	
end

%% Close database
SQL_closedatabase(dbc)



end