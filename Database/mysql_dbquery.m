function [outputData, errmessage] = mysql_dbquery(dbconn, sqlquery)
% mysql_dbquery     Retrieves data from a database connection.
%
% Can fetch data from either a JDBC4 connection, or using a Matlab database
% object (depending on the input class).

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
    curs = exec(dbconn,sqlquery);
    curs = fetch(curs);
    errmessage = curs.Message;

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
%         for j = 1:cols
%             qrfields{j} = char(rsmd.getColumnName(j));
%             if (nargout == 0)
%                 fprintf('%s\t', char(rsmd.getColumnName(j)));
%             end
%         end
%         if (nargout == 0)
%             fprintf('\n');
%         end
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
