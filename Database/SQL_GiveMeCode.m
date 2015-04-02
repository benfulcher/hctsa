% ------------------------------------------------------------------------------
% SQL_GiveMeCode
% ------------------------------------------------------------------------------
% Returns a string containing code for evaluating an operation with a given op_id.
% Can be difficult to do this manually, especially when dealing with structured
% outputs.
% 
%---INPUT:
% the_op_id, the operation_id you want code for.
% 
%---OUTPUT:
% The code as a function handle.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function TheCode = SQL_GiveMeCode(the_op_id);

% ------------------------------------------------------------------------------
% Open connection to database
[dbc, dbname] = SQL_opendatabase;

% Get opCode:
selectString = sprintf('SELECT Code FROM Operations WHERE op_id = %u',the_op_id);
opCode = mysql_dbquery(dbc,selectString);
opCode = opCode{1};

% Get MasterCode:
selectString = sprintf(['SELECT MasterLabel,MasterCode FROM MasterOperations WHERE mop_id = ' ...
                        '(SELECT mop_id FROM Operations WHERE op_id = %u)'],the_op_id);
mopData = mysql_dbquery(dbc,selectString);
mopLabel = mopData{1};
mopCode = mopData{2};

% Close connection to the mySQL database
SQL_closedatabase(dbc);
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Combine to give the required code as two function handles:
DotHere = regexp(opCode,'\.');
if ~isempty(DotHere) % A structure output
    WhatField = opCode(DotHere+1:end);
    % TheCode{1} = sprintf('@(x,y) %s;',mopCode); % Evaluate the structure
    TheCode = eval(sprintf('@(x,y) BF_GiveMeField(%s,''%s'');',mopCode,WhatField)); % Take the field from the structure
else
    TheCode = eval(sprintf('@(x,y) %s;',mopCode));
end
% ------------------------------------------------------------------------------

end