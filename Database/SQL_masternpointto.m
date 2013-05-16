%%% SQL_masternpointto
function SQL_masternpointto(dbname)
% Ben Fulcher 18/1/10 added dbname argument

if nargin<1
	dbname = '';
end

%% Open database
dbc = SQL_opendatabase(dbname);

	
%% Fills the NPointTo column in MasterOperations
selectstring = 'SELECT Master_id FROM MasterOperations';
[M_ids,qrf,rs,emsg] = mysql_dbquery(dbc,selectstring);
% Master ids, M_ids


if isempty(emsg)
	M_ids = vertcat(M_ids{:});
else
	disp('Error retrieving Master_ids');
end


for k=1:length(M_ids)
	updatestring = ['UPDATE MasterOperations SET NPointTo = ' ...
							'(SELECT COUNT(Master_id) FROM MasterPointerRelate WHERE Master_id = ' num2str(M_ids(k)) ') ' ...
					'WHERE Master_id = ' num2str(M_ids(k))];
	[rs,emsg] = mysql_dbexecute(dbc, updatestring);
	if ~isempty(emsg)
		disp(['Error with ' num2str(M_ids(k))]);
	end
end

end