function [qrcells, qrfields, queryresult, errmessage] = mysql_dbquery(dbconn, sqlquery)

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
    qrcells = {};
    i = 0;
    while (queryresult.next())
        i = i + 1;
        for j = 1:cols
            qrcells{i,j} = queryresult.getObject(j);
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
    qrfields = {};
    qrcells = {};
    queryresult = [];
    % le = lasterror;
    errmessage = le.message;
end

end