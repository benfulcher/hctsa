function [datafiles] = linsert(filename, handles,datafiles,s,nos)

% tstool/linsert
% status = linsert(filename, handles)
% Einfuegen eines Namens in die Liste, gleichzeitig wird die
% Selektion auf diesen Namen gesetzt

if isempty(s)
  [status,datafiles]=writelevel(filename,datafiles);
  [path,fname,ext,ver] = fileparts(filename);
else
  [status,datafiles]=siglevel(s,datafiles);
  fname=name(s);
end

if ~nos
  datafiles=sortdatafiles(datafiles);
  filllistbox(datafiles,handles,fname);
end