% Ben Fulcher 26/6/2013
% Move operations/masterops to new format where every operation links to a master operation

[dbc,dbname] = SQL_opendatabase;

% 1. Find all single operations in Operations Table
SelectString = 'SELECT OpName, Code, Keywords FROM Operations WHERE MasterLabel IS NULL';
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
sopname = qrc(:,1);
sopcode = qrc(:,2);
sopkw = qrc(:,3);
nsop = length(sopname);

% Check none have '.' in them
havedot = zeros(nsop,1);
for i = 1:nsop
    [a,b] = strtok(sopname{i},'.');
    if ~isempty(b)
        havedot(i) = 1;
    end
end
if sum(havedot) > 0
    input(sprintf('We found %u single operations with ''.'' in their name (show em?)...',sum(havedot)))
end
for i = 1:nsop
    if havedot(i)
        before = sopname{i};
        after = regexprep(before,'\.','');
        fprintf(1,'Reform %s to %s?\n',before,after);
        sopname{i} = after;
    end
end
input('What do you think?? (Cancel now or it will be too late)')


% 2. Find all pointers in Operations Table
SelectString = 'SELECT OpName, Code, Keywords FROM Operations WHERE MasterLabel IS NOT NULL';
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
popname = qrc(:,1);
popcode = qrc(:,2);
popkw = qrc(:,3);
npop = length(popname);

% 3. Retrieve all Master Operations
SelectString = 'SELECT MasterLabel, MasterCode FROM MasterOperations';
[qrc,~,~,emsg] = mysql_dbquery(dbc,SelectString);
moplab = qrc(:,1);
mopcode = qrc(:,2);
nmop = length(moplab);

%% Now, write back to input files
% 1. INP_ops_new
fid = fopen('Database/INP_ops_new.txt','w');
for i = 1:nsop
    fprintf(fid,'%s\t%s\t%s\n',sopname{i},sopname{i},sopkw{i});
end
for i = 1:npop
    fprintf(fid,'%s\t%s\t%s\n',popcode{i},popname{i},popkw{i});
end
fclose(fid);

% 2. INP_mops_new
fid = fopen('Database/INP_mops_new.txt','w');
for i = 1:nsop
    fprintf(fid,'%s\t%s\n',sopcode{i},sopname{i});
end
for i = 1:nmop
    fprintf(fid,'%s\t%s\n',mopcode{i},moplab{i});
end
fclose(fid);

fprintf(1,'Files written\n')

SQL_closedatabase(dbc);