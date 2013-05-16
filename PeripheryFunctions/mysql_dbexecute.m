function [execresult, errmessage] = mysql_dbexecute(dbconn, sqlcommand)
try
    dbstmt = dbconn.createStatement();
    execresult = dbstmt.execute(sqlcommand);
    errmessage = [];
catch
    le = lasterror;
    errmessage = le.message;
    execresult = [];
end
