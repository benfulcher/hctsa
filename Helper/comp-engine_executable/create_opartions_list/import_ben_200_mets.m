% ------------------------------------------------------------------------------
% import_ben_200_mets.m
% ------------------------------------------------------------------------------
% 
% Script to import The use 200 operations used for the pre-calculated
% feature vectors sent to me by ben in Oct 2014
% ------------------------------------------------------------------------------

% Reading in files for top 200 operations in old names and the map to the
% new names.

% call by fv_mets = import_ben_200_mets('TS_loc_guides_N.mat','oldnew_names_roslyn.csv');

% TODO implement Master function procedure. Now master functions are
% evaluated multiple times which creates processing overhead
function fv_mets = import_ben_200_mets(opNamesOldBen,oldnewnames)

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
%! create a cell array of the correct size to store the function calls
    fv_mets = cell(size(mlinkn));
%! create cell array for funtions whose call signature changed not only the
%  function name
    man_map = crt_old_new_mets_man();
    for i=1:size(mlinkn,1)
        %! this is a stand alone function and not dependent on a master
        %  operation
        if mlinkn(i) == 0
            % check if manual mapping is necessary
            ind = find(strcmpi(mcoden(i),man_map(:,1)));
            if ~isempty(ind)
                fv_mets(i) = man_map(ind,2);

            %! otherwise use automatic function name mapping and only change the 
            % function name
            else
                replaced_str = repl_comm(mcoden(i),name_map_in);
                if isempty(replaced_str)
                    fv_mets(i) = 'NaN';
                else 
                    fv_mets(i) = replaced_str;
                end
            end
        else
            tmp = mcoden(i);
            [~,pointer] = strtok(mcoden(i),'.');
            % check if manual mapping is necessary
            ind = find(strcmpi(strcat(Mmcode(mlinkn(i)),pointer),man_map(:,1)));
            if ~isempty(ind)
                fv_mets(i) = man_map(ind,2);  
            else
                % otherwise use automatic function name mapping            
                replaced_str = repl_comm(Mmcode(mlinkn(i)),name_map_in);
                if isempty(replaced_str)
                    fv_mets(i) = 'NaN';
                else 
                    fv_mets(i) = strcat(replaced_str,pointer);
                end
            end            
        end
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
