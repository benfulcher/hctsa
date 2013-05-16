function makedir(dirname)

% create directory dirname

[path,name] = fileparts(dirname);

if ~exist(path,'dir') 
    makedir(path);
end

if isunix
    if unix(['mkdir ' dirname])~=0
        error(['Cannot create directory ' dirname])
    end
else
    mkdir(path, name);
end
