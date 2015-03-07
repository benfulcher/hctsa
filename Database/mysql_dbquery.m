% ------------------------------------------------------------------------------
% mysql_dbquery
% ------------------------------------------------------------------------------
% Used to retrieve data from a database connection
% ------------------------------------------------------------------------------
%---HISTORY:
% Ben Fulcher, 2015-01-30. Adds support to fetch data using a Matlab database object.
% Removed the additional outputs that are never actually used.
% ------------------------------------------------------------------------------

function [outputData, errmessage] = mysql_dbquery(dbconn, sqlquery)

switch class(dbconn)
case 'database'
    % Database connection opened using database toolbox:
    curs = exec(dbconn,sqlquery);
    curs = fetch(curs);
    errmessage = curs.Message;
    qrfields = {};
    queryresult = [];
    
    if isempty(curs.Message)
        % Success!:
        outputData = curs.Data;
    else
        outputData = {};
    end
    
    % Close cursor object:
    close(curs);
    
    % If 'No Data', change to empty, to match the syntax of 'com.mysql.jdbc.JDBC4Connection'
    if isempty(errmessage) && iscell(outputData) && ischar(outputData{1}) && strcmp(outputData{1},'No Data')
        outputData = {};
    end
    
case 'com.mysql.jdbc.JDBC4Connection'
    % Database connection opened using java commands:
    try
        dbstmt = dbconn.createStatement();
        queryresult = dbstmt.executeQuery(sqlquery);
        rsmd = queryresult.getMetaData();
        cols = rsmd.getColumnCount();
        for j = 1:cols
            qrfields{j} = char(rsmd.getColumnName(j));
            if (nargout == 0)
                fprintf('%s\t', char(rsmd.getColumnName(j)));
            end
        end
        if (nargout == 0)
            fprintf('\n');
        end
        outputData = {};
        i = 0;
        while (queryresult.next())
            i = i + 1;
            for j = 1:cols
                outputData{i,j} = queryresult.getObject(j);
                if (nargout == 0)
                    fprintf('%s\t', char(queryresult.getString(j)));
                end
            end
            if (nargout == 0)
                fprintf('\n');
            end
        end
        errmessage = [];
    catch le
        % qrfields = {};
        outputData = {};
        % queryresult = [];
        % le = lasterror;
        errmessage = le.message;
    end
    
otherwise
    error('Unknown database connection class %s',class(dbc));
end


end