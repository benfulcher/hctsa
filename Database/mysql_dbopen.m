function [dbconnection, errmessage] = mysql_dbopen(serverhost, dbname, uname, pword)

import java.lang.Thread;
import java.lang.Class;
import java.sql.DriverManager;

ct = java.lang.Thread.currentThread();
cl = ct.getContextClassLoader();

try
    java.lang.Class.forName('com.mysql.jdbc.Driver', true, cl);
    dburl = sprintf('jdbc:mysql://%s/%s', serverhost, dbname);
    dbconnection = java.sql.DriverManager.getConnection(dburl, uname, pword);
    errmessage = [];
catch le
%     le = lasterror;
    errmessage = le.message;
    disp(errmessage);
    dbconnection = [];
end
