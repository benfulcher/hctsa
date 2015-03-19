ind_arr = [1:20];

src_folder = '/home/philip/work/Git_repositories/hctsa/Helper/speciesID/Datarfb/hctsa_final/';
cd /home/philip/work/Git_repositories/hctsa/Helper/speciesID/Datarfb/hctsa_final/;

is_first = true;
for i=ind_arr
    i
    filename = [src_folder sprintf('HCTSA_loc_%03d.mat',i)];
    if is_first
        copyfile(filename, [src_folder 'HCTSA_loc_combined.mat']);
        is_first = false;
    else
        TSQ_combine('HCTSA_loc_combined.mat',filename)
    end
end
