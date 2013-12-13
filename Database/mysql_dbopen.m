function [dbconnection, errmsg] = mysql_dbopen(serverhost, dbname, uname, pword)

import java.lang.Thread;
import java.lang.Class;
import java.sql.DriverManager;

ct = java.lang.Thread.currentThread();
cl = ct.getContextClassLoader();

% Open database
errmsg = []; % error message

% First check driver
try
    java.lang.Class.forName('com.mysql.jdbc.Driver', true, cl);
catch le
    errmsg = le.message;
    dbconnection = [];
    error('Error with java database connector: %s',errmsg);
    % actually return an error
end

% Now try to connect
try
    dburl = sprintf('jdbc:mysql://%s/%s', serverhost, dbname);
    dbconnection = java.sql.DriverManager.getConnection(dburl, uname, pword);
catch le
    fprintf(1,'Error connecting to the database ''%s'' at ''%s'':\n%s\n',dbname,serverhost,le.message);
    fprintf(1,['Perhaps an incorrect username (''%s'') and password (''%s'') combination?\n'], uname, pword);
    dbconnection = [];
    % Not really a Matlab 'error' -- just print the suspected problem to screen
end

end