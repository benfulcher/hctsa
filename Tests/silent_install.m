% INSTALL   Installs the hctsa code package from scratch.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------
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

%-------------------------------------------------------------------------------
%
fprintf('All done!')
