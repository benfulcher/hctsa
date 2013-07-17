function SQL_linkpointermaster(dbname)
% Ben Fulcher 18/1/10 Added dbname argument

if nargin < 1
	dbname = '';
end

%% Open database
[dbc,dbname] = SQL_opendatabase(dbname);

%% Creates the MasterPointerRelate Table by linking metrics with their masters
% Look through code column of operations and retrieve bits before "."
% These must be the master structure that should link
SelectString = 'SELECT m_id, Code FROM Operations WHERE Pointer = 1';
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,SelectString);
pointermids = vertcat(qrc{:,1}); % column vector
pointercodes = qrc(:,2); % column cell of strings
numpointers = length(pointercodes); % number of pointer metrics that need a link
pointerrefs = cell(numpointers,1);

% Filter off the bits after '.' -- store in pointerrefs
for i = 1:numpointers
	dotishere = strfind(pointercodes{i},'.');
	pointerrefs{i} = pointercodes{i}(1:dotishere-1);
end

% Create new linking table if it doesn't already exist (in which case, delete it and start again)
[qrc,qrf,rs,emsg] = mysql_dbquery(dbc,'Show Tables');

if ismember('MasterPointerRelate',qrc)
	% already exists -- drop the table then recreate it
	[rs,emsg] = mysql_dbexecute(dbc, 'DROP TABLE MasterPointerRelate');
	if ~isempty(emsg)
		disp(['Error dropping MasterPointerRelate table from ' dbname]);
		disp(emsg); beep
		keyboard
	end
end
createstring = ['CREATE TABLE MasterPointerRelate (mop_id integer, m_id integer, ' ...
				'FOREIGN KEY (m_id) REFERENCES Operations(m_id) ON DELETE CASCADE ON UPDATE CASCADE, ' ...
				'FOREIGN KEY (mop_id) REFERENCES MasterOperations(mop_id) ON DELETE CASCADE ON UPDATE CASCADE)'];
[rs,emsg] = mysql_dbexecute(dbc, createstring);
if isempty(emsg)
	disp(['Created MasterPointerRelate table in ' dbname]);
else
	disp('Error creating table: MasterPointerRelate');
	disp(emsg); beep
	keyboard
end


% Write relationships
disp(['Writing relationships between master and pointer operations to MasterPointerRelate table in database. ' ...
			'\n Be patient, this may take some minutes...']);
times = zeros(numpointers,1);
for i = 1:numpointers
	tic
	insertstring = ['INSERT INTO MasterPointerRelate (m_id, mop_id) SELECT ' num2str(pointermids(i)) ', mop_id FROM MasterOperations' ...
	    	' WHERE (MasterLabel = ''' pointerrefs{i} ''')'];
	[rs,emsg] = mysql_dbexecute(dbc, insertstring);
	if ~isempty(emsg)
		disp(['Error executing MasterPointer relationship with ' pointerrefs{i}]);
		disp([insertstring]); disp(emsg); beep
		keyboard
	end
	
	times(i) = toc;
	if mod(i,floor(numpointers/10))==0
		disp(['Approximately ' BF_thetime(mean(times(1:i))*(numpointers-i)) ' remaining!']);
	end
end
end