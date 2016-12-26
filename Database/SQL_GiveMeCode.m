function theCode = SQL_GiveMeCode(op_id)
% SQL_GiveMeCode a string containing code for evaluating an operation with a given op_id
%
% Can be difficult to do this manually, especially when dealing with structured
% outputs.
%
%---INPUT:
% op_id, the operation_id you want code for.
%
%---OUTPUT:
% The code as a function handle.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Open connection to database
[dbc, ~] = SQL_opendatabase;

% Get opCode:
selectString = sprintf('SELECT Code FROM Operations WHERE op_id = %u',op_id);
opCode = mysql_dbquery(dbc,selectString);
opCode = opCode{1};

% Get MasterCode:
selectString = sprintf(['SELECT MasterCode FROM MasterOperations WHERE mop_id = ' ...
                        '(SELECT mop_id FROM Operations WHERE op_id = %u)'],op_id);
mopData = mysql_dbquery(dbc,selectString);
mopCode = mopData{1};

% Close connection to the mySQL database
SQL_closedatabase(dbc);
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Combine to give the required code as two function handles:
dotHere = regexp(opCode,'\.');
if ~isempty(dotHere) % A structure output
    whatField = opCode(dotHere+1:end);
    theCode = eval(sprintf('@(x,y) BF_GiveMeField(%s,''%s'');',mopCode,whatField)); % Take the field from the structure
else
    theCode = eval(sprintf('@(x,y) %s;',mopCode));
end
% ------------------------------------------------------------------------------

end
