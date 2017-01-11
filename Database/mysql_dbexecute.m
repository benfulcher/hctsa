function [execResult, errMessage] = mysql_dbexecute(dbconn, sqlCommand)
% MYSQL_DBEXECUTE
%
% Retrieves data from a database connection, using either the Matlab
% database toolbox, or a java JDBC4 connection (depending on the input
% connection class).

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

switch class(dbconn)
case 'database'
    % Database connection opened using database toolbox:
    curs = exec(dbconn,sqlCommand);
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
        execResult = dbstmt.execute(sqlCommand);
        errMessage = [];
    catch emsg
        errMessage = emsg.message;
        execResult = [];
    end

otherwise
    error('Unknown database connection class %s',class(dbconn));
end

end
