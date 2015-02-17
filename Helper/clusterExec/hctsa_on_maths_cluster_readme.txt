1.) Create as many different HCTSA_loc.mat files as jobs you want to run in parrallel (containing ts and op you want to run in the respective job)
2.) Put those HCTSA_loc.mat files on the cluster (e.g. /scratchcomp10/pknaute/HCTSA_loc_1/HCTSA_loc_001.mat ... )
3.) Create a plain-text file called 'file_list.text' containing the full paths of all HCTSA_loc.mat files you want to calculate. Each path a single line in the file terminated by a 'newline'
4.) Put the file_list.text,submit_array.pbs and simple_calc.m in the same folder on the cluster (doesn't have to be the same folder as the HCTSA_loc.mat files)
5.) If you want to run e.g. 10 (serial) jobs in paralell run 
'qsub -t 1-10 submit_array.pbs' while in the folder containing the .pbs file. To be able to submit jobs you have to be logged in on a macomp01 .. macomp16, the macomp000 won't let you (at least not me).

