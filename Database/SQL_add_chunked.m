% ------------------------------------------------------------------------------
% SQL_add_chunked
% ------------------------------------------------------------------------------
% 
% Insert a set of things into the database
% insertstring is the insert portion of the query
% dataset is a cell array of formatted strings like {'(''abc'',1)'}
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Romesh Abeysuriya, February 2013
% Ben D. Fulcher <ben.d.fulcher@gmail.com>, <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function SQL_add_chunked(dbc,insertstring,dataset,isduplicate,chunksize)

if nargin < 4 || isempty(isduplicate)
    isduplicate = zeros(size(dataset));
end

if nargin < 5 || isempty(chunksize)
    chunksize = 100; % run this many queries at a time
    % This parameter can be tweaked depend on the value of max_allowed_packet
    % on the mySQL server.
end

for k = 1:chunksize:length(dataset)
    query = insertstring; % start with the insert statement
    for j = k:min(k+chunksize-1,length(dataset)) % don't duplicate
        if ~isduplicate(j)
            query = sprintf('%s %s,',query,dataset{j}); % add values in parentheses in dataset{j}
        end
    end
    query = query(1:end-1); % remove the final comma
    
    [~, emsg] = mysql_dbexecute(dbc,query); % evaluate this chunk
    if ~isempty(emsg)
        fprintf(1,'Error in SQL_add_chunked...\n%s\n',emsg)
        keyboard
    end
end

end