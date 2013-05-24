function emsg = SQL_closedatabase(dbc)
% Closes the connection to database dbc

try
   dbc.close();
catch emsg
	disp('Error closing database')
end

clear dbc

end