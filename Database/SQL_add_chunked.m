function SQL_add_chunked(dbc,insertstring,dataset,isduplicate,append)
    % Insert a set of things into the database
    % insertstring is the insert portion of the query
    % dataset is a cell array of formatted strings like {'(''abc'',1)'}
    % Romesh Abeysuriya, February 2013
    
    if nargin < 4 || isempty(isduplicate)
        isduplicate = zeros(size(dataset));
    end
    
    if nargin < 5 || isempty(append)
        append = '';
    end
        
    chunksize = 1000;

    for k = 1:chunksize:length(dataset)      
        query = insertstring;
        for j = k:k+chunksize
            if j > length(dataset)
                break
            elseif ~isduplicate(j)
                query = sprintf('%s %s,',query,dataset{j});
            end
        end
        [rs,emsg] = mysql_dbexecute(dbc,[query(1:end-1) ' ' append]);
        if ~isempty(emsg)
            fprintf('Error in SQL_add_chunked...?!')
            keyboard
        end
    end
    
end