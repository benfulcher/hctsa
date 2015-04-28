% ------------------------------------------------------------------------------
% SQL_add_chunked
% ------------------------------------------------------------------------------
% 
% Insert a large set of time series or operations into the database using
% repeated queries, adding smaller subsets over multiple iterations.
% 
%---INPUTS:
% 
% dbc, the database connection
% 
% insertString, the insert portion of the query
% 
% dataset, a cell array of formatted strings like {'(''abc'',1)'}
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013, Romesh Abeysuriya
% Ben D. Fulcher <ben.d.fulcher@gmail.com>, <http://www.benfulcher.com>
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

function SQL_add_chunked(dbc,insertString,dataset,isDuplicate,chunkSize)

% ------------------------------------------------------------------------------
% Check inputs
% ------------------------------------------------------------------------------
if nargin < 4 || isempty(isDuplicate)
    isDuplicate = zeros(size(dataset));
end

if nargin < 5 || isempty(chunkSize)
    chunkSize = 50; % Run this many queries at a time
    % This parameter can be tweaked depend on the value of max_allowed_packet
    % on the mySQL server.
end

% ------------------------------------------------------------------------------
% Start adding chunks to the database
% ------------------------------------------------------------------------------
for k = 1:chunkSize:length(dataset)
    query = insertString; % Start with the insert statement
    for j = k:min(k+chunkSize-1,length(dataset)) % Don't repeat statements
        if ~isDuplicate(j)
            query = sprintf('%s %s,',query,dataset{j}); % Add values in parentheses in dataset{j}
        end
    end
    
    if (length(query) > length(insertString))
        % There's something to be added, i.e., not all isDuplicate in this chunk:
    
        query = query(1:end-1); % Remove the final comma
    
        [~, emsg] = mysql_dbexecute(dbc,query); % Evaluate this chunk
        if ~isempty(emsg)
            error('Error in SQL_add_chunked for chunk %u with chunk size %u...\n%s\n',k,chunkSize,emsg)
        end
    end
end

end