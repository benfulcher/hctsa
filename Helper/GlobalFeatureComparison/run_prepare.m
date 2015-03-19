
% ------------------------------------------------------------------------------
% run_prepare
% ------------------------------------------------------------------------------
%
% prepares a HCTSA_loc.mat file for external timeseries (not included in
% the database) using the PK_prepare_external_ts() function.
% 
% ------------------------------------------------------------------------------

%opids = SQL_getids('ops',1,{},{'shit','kalafutvisscher','waveletTB'});
%PK_prepare_external_ts(1);
ts_id = 1;
root_dir = '/home/philip/work/OperationImportanceProject/dataset/';
sub_dirs = dir(root_dir);
sub_dirs = sub_dirs(3:end);
for i=1:size(sub_dirs)
    sub_dir_name = [root_dir '/' sub_dirs(i).name];
    files = dir(sub_dir_name);
    files = files(3:end);
    for k = 1:size(files)
        file_name = [sub_dir_name '/' files(k).name];
        file_data = load(file_name);
        for l=1:size(file_data,1)
            if ts_id == 1
                ts_info = {ts_id,file_name,[int2str(l) ',' int2str(file_data(l,1)) ',' ], ...
                            size(file_data,2)-1,file_data(l,2:end) };
                ts_id = ts_id + 1;
            else
                ts_info = [ts_info; {ts_id,file_name,[int2str(l) ',' int2str(file_data(l,1)) ',' ], ...
                           size(file_data,2)-1,file_data(l,2:end) } ];
                ts_id = ts_id + 1;
            end
        end
    end
    
end
