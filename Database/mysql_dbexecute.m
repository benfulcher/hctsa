% ------------------------------------------------------------------------------
% mysql_dbquery
% ------------------------------------------------------------------------------
% Used to retrieve data from a database connection
% ------------------------------------------------------------------------------
%---HISTORY:
% Ben Fulcher, 2015-01-30. Adds support to fetch data using a Matlab database object.
% Removed the additional outputs that are never actually used.
% ------------------------------------------------------------------------------

function [execResult, errMessage] = mysql_dbexecute(dbconn, sqlcommand)

switch class(dbconn)
case 'database'
    % Database connection opened using database toolbox:
    curs = exec(dbconn,sqlcommand);
    errMessage = curs.Message;

    % Close cursor object:
    close(curs);
    
    if isempty(errMessage)
        % Success!:
        execResult = 0;
    else
        execResult = [];
    end

case 'com.mysql.jdbc.JDBC4Connection'
    % Database connection opened using java commands:
    try
        dbstmt = dbconn.createStatement();
        execResult = dbstmt.execute(sqlcommand);
        errMessage = [];
    catch
        le = lasterror;
        errMessage = le.message;
        execResult = [];
    end

otherwise
    error('Unknown database connection class %s',class(dbc));
end

end