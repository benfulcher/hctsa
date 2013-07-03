function SQL_add_chunked(dbc,insertstring,dataset,isduplicate)
% Insert a set of things into the database
% insertstring is the insert portion of the query
% dataset is a cell array of formatted strings like {'(''abc'',1)'}
% Romesh Abeysuriya, February 2013

if nargin < 4 || isempty(isduplicate)
    isduplicate = zeros(size(dataset));
end
    
chunksize = 1000;

for k = 1:chunksize:length(dataset)
    query = insertstring; % start with the insert statement
    for j = k:min(k+chunksize-1,length(dataset)) % don't duplicate
        if ~isduplicate(j)
            query = sprintf('%s %s,',query,dataset{j}); % add values in parentheses in dataset{j}
        end
    end
    query = query(1:end-1); % remove the final comma
    
    [~,emsg] = mysql_dbexecute(dbc,query); % evaluate this chunk
    if ~isempty(emsg)
        fprintf(1,'Error in SQL_add_chunked...\n%s\n',emsg)
        keyboard
    end
end

end