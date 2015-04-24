% ------------------------------------------------------------------------------
% mysql_dbopen
% ------------------------------------------------------------------------------
% Opens a connection to the database using the mySQL j-connector.
% Checks for an available database toolbox and uses that, but otherwise uses
% java commands that should achieve the same thing.
% 
%---HISTORY:
% Ben Fulcher, 2014-12-12: Some users reported issues using this function,
% changed to use the database toolbox by default (if a licence exists)
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

function [dbConnection, errMessage] = mysql_dbopen(serverHost, databaseName, userName, password, customPort, useDBToolbox)

% Set defaults:
if nargin < 5
    customPort = 3306; % use default port for mySQL: 3306
end
if nargin < 6
    useDBToolbox = 0; % This is much slower than using java directly
end

% ------------------------------------------------------------------------------
% Connect to the database (using either the database toolbox or the mySQL J-connector)
% ------------------------------------------------------------------------------

% Check if a Matlab database toolbox exists, and use it if specified:
dbToolboxExists = license('test','database_toolbox');

if useDBToolbox && dbToolboxExists
    % ------------------------------------------------------------------------------
    % Connect using the Matlab's database toolbox
    % ------------------------------------------------------------------------------
    
    % The JDBC interface is slower but platform independent
    dbConnection = database(databaseName,userName,password,'Vendor','MySQL','Server',serverHost,'PortNumber',customPort);

    errMessage = dbConnection.Message;
    
    % Check for a connection error:
    if ~isempty(errMessage)
        error(['Error connecting to the database using Matlab database toolbox:\n' ...
                        '%s (pass:%s) connecting to %s at %s\n%s\n'],  ...
                        userName,password,databaseName,serverHost,dbConnection.Message);
    end
    
else
    % ------------------------------------------------------------------------------
    % Connect using the mySQL J-connector (if no database toolbox exists)
    % ------------------------------------------------------------------------------
    import java.lang.Thread;
    import java.lang.Class;
    import java.sql.DriverManager;

    ct = java.lang.Thread.currentThread();
    cl = ct.getContextClassLoader();
    % ------------------------------------------------------------------------------

    % First check driver:
    try
        java.lang.Class.forName('com.mysql.jdbc.Driver', true, cl);
    catch le
        errMessage = le.message;
        dbConnection = [];
        error('Error with java database connector:\n%s',errMessage);
    end

    % ------------------------------------------------------------------------------
    % Now try to connect:
    try
        dburl = sprintf('jdbc:mysql://%s/%s', serverHost, databaseName);
        dbConnection = java.sql.DriverManager.getConnection(dburl, userName, password);
    catch le
        dbConnection = [];
        error('Error connecting to the database ''%s'' at ''%s'':\n%s\n',databaseName,serverHost,le.message);
        % fprintf(1,['\nPerhaps due to an incorrect username (''%s'') and password (''%s'') combination?\n'], userName, password);
    end
end

end
