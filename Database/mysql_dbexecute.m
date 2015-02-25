% ------------------------------------------------------------------------------
% mysql_dbquery
% ------------------------------------------------------------------------------
% Used to retrieve data from a database connection
% ------------------------------------------------------------------------------
%---HISTORY:
% Ben Fulcher, 2015-01-30. Adds support to fetch data using a Matlab database object.
% Removed the additional outputs that are never actually used.
% ------------------------------------------------------------------------------

function [execresult, errmessage] = mysql_dbexecute(dbconn, sqlcommand)

switch class(dbconn)
case 'database'
    % Database connection opened using database toolbox:
    curs = exec(dbconn,sqlcommand);
    errmessage = curs.Message;

    % Close cursor object:
    close(curs);
    
    if isempty(errmessage)
        % Success!:
        execresult = 0;
    else
        execresult = [];
    end

case 'com.mysql.jdbc.JDBC4Connection'
    % Database connection opened using java commands:
    try
        dbstmt = dbconn.createStatement();
        execresult = dbstmt.execute(sqlcommand);
        errmessage = [];
    catch
        le = lasterror;
        errmessage = le.message;
        execresult = [];
    end

otherwise
    error('Unknown database connection class %s',class(dbc));
end

end