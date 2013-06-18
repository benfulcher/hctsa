function SQL_masternpointto(dbname)
% Populates the NPointTo column in MasterOperations
% NPointTo gives the number of 
% Ben Fulcher 18/1/10 added dbname argument

if nargin < 1
	dbname = ''; % use the default database
end

%% Open database
dbc = SQL_opendatabase(dbname);
	
%% Fills the NPointTo column in MasterOperations
SelectString = 'SELECT mop_id FROM MasterOperations';
[M_ids,~,~,emsg] = mysql_dbquery(dbc,SelectString);

if isempty(emsg)
	M_ids = vertcat(M_ids{:}); % vector of master_ids
else
	error('SQL_masternpointto: Error retrieving mop_ids');
end


for k = 1:length(M_ids)
	UpdateString = ['UPDATE MasterOperations SET NPointTo = ' ...
							'(SELECT COUNT(mop_id) FROM MasterPointerRelate WHERE mop_id = ' num2str(M_ids(k)) ') ' ...
					'WHERE mop_id = ' num2str(M_ids(k))];
	[rs,emsg] = mysql_dbexecute(dbc, UpdateString);
	if ~isempty(emsg)
		disp(['Error with ' num2str(M_ids(k))]);
	end
end

end