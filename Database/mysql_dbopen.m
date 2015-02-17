function [dbconnection, errmsg] = mysql_dbopen(serverhost, dbname, uname, pword)

import java.lang.Thread;
import java.lang.Class;
import java.sql.DriverManager;

ct = java.lang.Thread.currentThread();
cl = ct.getContextClassLoader();

% ------------------------------------------------------------------------------

errmsg = ''; % error message empty by default

% % First check driver:
% try
%     java.lang.Class.forName('com.mysql.jdbc.Driver', true, cl);
% catch le
%     errmsg = le.message;
%     dbconnection = [];
%     error('Error with java database connector: %s',errmsg);
% end

% ------------------------------------------------------------------------------
% Now try to connect:
try
    dburl = sprintf('jdbc:mysql://%s/%s', serverhost, dbname);
    dbconnection = java.sql.DriverManager.getConnection(dburl, uname, pword);
catch le
    dbconnection = [];
    error('Error connecting to the database ''%s'' at ''%s'':\n%s\n',dbname,serverhost,le.message);
    % fprintf(1,['\nPerhaps due to an incorrect username (''%s'') and password (''%s'') combination?\n'], uname, pword);
end

end