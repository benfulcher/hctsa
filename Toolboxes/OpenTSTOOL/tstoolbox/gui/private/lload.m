function [datafiles] = lload(filename, handles,datafiles)

% tstool/lload


[status,datafiles]=writelevel(filename,datafiles);
[path,fname,ext,ver] = fileparts(filename);


datafiles=sortdatafiles(datafiles)
filllistbox(datafiles,handles,fname);
