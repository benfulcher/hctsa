%%% SQL_closedatabase
function errmessage = SQL_closedatabase(dbc)
% Closes the connection to database dbc

try
   dbc.close();
catch emsg
	disp('Error closing database')
end

% if isempty(emsg)
% 	disp('Database closed')
% end

clear dbc

end