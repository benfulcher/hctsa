function SQL_closedatabase(dbc)
% SQL_closedatabase     Closes the connection to database, dbc.

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

switch class(dbc)
case 'com.mysql.jdbc.JDBC4Connection'
    % Database connection opened using java commands:
    dbc.close();
case 'database'
    % Database connection opened using database toolbox:
    close(dbc);
otherwise
    error('Unknown database connection class %s',class(dbc));
end


end
