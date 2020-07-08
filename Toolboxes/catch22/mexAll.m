% path where C implementation lives
basePath = './C';
ipath = {['-I', basePath], '-I.'};

% list all C files to include in the mex-call
CFiles = dir(basePath);
CFiles = CFiles(cellfun(@(x) contains(x, '.c'), {CFiles.name}));
CFileNames = {CFiles.name};

% add path to the c file names
includeFiles = cellfun(@(x) fullfile(basePath, x), CFileNames, 'UniformOutput', false);

% get function names
fid = fopen('featureList.txt','r');
i = 1;
tline = fgetl(fid);
featureNames{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    featureNames{i} = tline;
end
fclose(fid);

% mex all feature functions separately
for i = 1:length(featureNames)-1

    featureName = featureNames{i};

    fprintf('Compiling %s...\n', featureName);
    mex(ipath{:}, ['catch22_', featureName,'.c'], 'M_wrapper.c', includeFiles{:})
    fprintf('\n');

end
