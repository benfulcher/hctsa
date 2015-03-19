% Remember home directory
MyPBSHome = pwd

% Load paths for the HCTSA package
cd ~/work/matlab/hctsa/
startup

% Move Matlab back to the working PBS directory
cd(MyPBSHome)

% check if a custom TSQ_loc filename has been defined.
if not (exist('HCTSA_loc_filename'))
    HCTSA_loc_filename = 'HCTSA_loc.mat'
end

TSQ_brawn(0,0,1,HCTSA_loc_filename); 

clear
exit
