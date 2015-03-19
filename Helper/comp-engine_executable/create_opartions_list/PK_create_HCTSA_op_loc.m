function op_id = PK_create_HCTSA_op_loc(opNamesOldBen,oldnewnames)

%! load the operation commands and names from the precalculated
%  TS_loc_guides_N.mat which was provided by ben
    ben_200_op_old = load(opNamesOldBen);

    mcoden = ben_200_op_old.mcoden;
    mlinkn = ben_200_op_old.mlinkn;
    Mmcode = ben_200_op_old.Mmcode;
%! load the map from old to new names which is given in a (incomplete) .csv file
%  by Roslyn.
    file_in_name = oldnewnames;
    fid = fopen(file_in_name);
    name_map_in = textscan(fid,'%s %s %s','delimiter',',','CollectOutput',1);
    name_map_in = name_map_in{1}; 
    name_map_in = name_map_in(2:end,2:end); 
    fclose(fid);
    
    [dbc, ~] = SQL_opendatabase();
    %! create an array of the correct size to store the op_ids
    op_id = zeros(200,1);
    %! create cell array for funtions whose call signature changed not only the
    %  function name
    man_map = crt_old_new_mets_man();
    for i=1:size(mlinkn,1)
        %! this is a stand alone function and not dependent on a master
        %  operation
        if mlinkn(i) == 0
            pointer = {''};
            % check if manual mapping is necessary
            ind = find(strcmpi(mcoden(i),man_map(:,1)));
            if ~isempty(ind)
                replaced_str = man_map(ind,2);
                pointer = man_map(ind,3);
            %! otherwise use automatic function name mapping and only change the 
            % function name
            else
                replaced_str = repl_comm(mcoden(i),name_map_in);
            end
        else
            [~,pointer] = strtok(mcoden(i),'.');
            % check if manual mapping is necessary
            ind = find(strcmpi(strcat(Mmcode(mlinkn(i)),pointer),man_map(:,1)));
            if ~isempty(ind)
                replaced_str = man_map(ind,2);
                pointer = man_map(ind,3);
                %replaced_str = man_map(ind,2);
                %pointer = man_map(ind,3);
            else
                % otherwise use automatic function name mapping 
                replaced_str = repl_comm(Mmcode(mlinkn(i)),name_map_in);
            end  
        end
        
        % --------------------------------------------------------------------
        op_command = strrep(replaced_str{1}, '''', '\'''); 
        sql_template = ['SELECT op_id FROM Operations AS Op JOIN MasterOperations AS MOp ON Op.mop_id = MOp.mop_id ' ...
                'WHERE MOp.MasterCode LIKE  ''%s''AND Op.Code LIKE  ''%%%s%%'' '];
        sqlcommand = sprintf(sql_template,op_command,pointer{1});
        execresult = mysql_dbquery(dbc, sqlcommand);
        % --------------------------------------------------------------------              
        
        
        % --------------------------------------------------------------------
        if ~(isempty(execresult))
            op_id(i) = execresult{1};
        end
        % --------------------------------------------------------------------
    end
end

%! replace the old function name with the new function name.
% \param old_comm [string] old command
% \param comm_map [cell array] mapping between old and new function names
% \returns [string] the function call with new function name.
function replaced_str = repl_comm(old_comm,comm_map)
    % split the old command in name and arguments
    [comm_l,comm_r] = strtok(old_comm,'(');
    % find the index of the old command in the commands map
    ind = find(strcmpi(comm_l,comm_map(:,1)));
    % combine new function name with parameter values to a new string
    replaced_str = strcat(comm_map(ind,2),comm_r);
end
    