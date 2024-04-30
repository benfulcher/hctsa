% modified version of the original install.m, designed to run from the
% command line in a Git Actions runner. 

cd('../') 
% ------------------------------------------------------------------------------
% 0. Check that the latest git submodules are loaded (catch22)
% ------------------------------------------------------------------------------
if exist(fullfile(pwd,'Toolboxes','catch22','wrap_Matlab'),'dir')==0
    system('git submodule update --init --recursive');
end

% ------------------------------------------------------------------------------
%% 1. Add the paths:
% ------------------------------------------------------------------------------
fprintf(1,'-1- Adding paths needed for the repository...\n');
try
	startup
catch emsg
	fprintf(1,'error.\n%s\n',emsg.message);
end
fprintf(1,'\n');

% ------------------------------------------------------------------------------
%% 2. Check the toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,'-2- Checking installation of relevant Matlab Toolboxes...\n');

% Essential:
noStatToolbox = BF_CheckToolbox('statistics_toolbox',true,true);
if noStatToolbox
    fprintf(1,'The Statistics and Machine Learning Toolbox is required for hctsa and needs to be installed.\n');
end
% For distributed computing
noDistributedToolbox = BF_CheckToolbox('parallel_computing_toolbox',true,true);
if noDistributedToolbox
    fprintf(1,'Matlab''s Parallel Computing Toolbox is required for distributed calculation but is not installed.\n');
end
% Desirable (needed for some hctsa features):
toolboxCodes = {'curve_fitting_toolbox','signal_toolbox','identification_toolbox',...
                    'wavelet_toolbox','econometrics_toolbox','financial_toolbox'};
flags = zeros(8,1);
names = cell(8,1);
for i = 1:length(toolboxCodes)
    [flags(i),names{i}] = BF_CheckToolbox(toolboxCodes{i},true,true);
end

fprintf(1,'%u/%u hctsa-relevant Matlab toolboxes are installed.\n',sum(~flags),length(flags));
if any(flags)
    fprintf(1,'For full functionality, please install:\n');
    theNo = find(flags);
    numNo = length(theNo);
    for i = 1:numNo
        fprintf(1,'%s\n',names{theNo(i)});
    end
else
    fprintf(1,'Looking good for Matlab toolbox installation! :)\n');
end
fprintf(1,'\n');

% ------------------------------------------------------------------------------
%% 3. Attempt to compile the executables required by the periphery Toolboxes:
% ------------------------------------------------------------------------------
cd Toolboxes
compile_mex
cd('../');
fprintf('All done!')
