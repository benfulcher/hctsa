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

fprintf(1,['This script will set up the Highly Comparative Time-Series ' ...
                        'Analysis (hctsa) code package from scratch!\n']);
fprintf(1,['We will:' ...
            '\n-1- Add the paths needed for the repository,' ...
            '\n-2- Check toolboxes,' ...
            '\n-3- Compile the external time-series toolboxes for this system.\n\n']);

% ------------------------------------------------------------------------------
% 0. Check that the latest git submodules are loaded (catch22)
% ------------------------------------------------------------------------------
if exist(fullfile(pwd,'Toolboxes','catch22','wrap_Matlab'),'dir')==0
    input('You may not have added the git submodule for catch22... Press any key to attempt this now...')
    system('git submodule update --init --recursive');
end

% ------------------------------------------------------------------------------
%% 1. Add the paths:
% ------------------------------------------------------------------------------
fprintf(1,'-1- Adding paths needed for the repository...\n');
input('<<Press any key to continue>>')
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
noStatToolbox = BF_CheckToolbox('statistics_toolbox',true);
if noStatToolbox
    fprintf(1,'The Statistics and Machine Learning Toolbox is required for hctsa and needs to be installed.\n');
end
% Desirable:
toolboxCodes = {'curve_fitting_toolbox','signal_toolbox','identification_toolbox',...
                    'wavelet_toolbox','econometrics_toolbox','robust_toolbox',...
                    'financial_toolbox','distrib_computing_toolbox'};
flags = zeros(8,1);
names = cell(8,1);
for i = 1:length(toolboxCodes)
    [flags(i),names{i}] = BF_CheckToolbox(toolboxCodes{i},true);
end

fprintf(1,'%u/%u hctsa-relevant Matlab toolboxes are installed.\n',sum(flags),length(flags));
if any(flags)
    fprintf(1,'For full functionality, please install:\n');
    theNo = find(flags);
    numNo = length(theNo);
    for i = 1:numNo
        fprintf(1,'%s\n',names{theNo(i)});
    end
else
    fprintf(1,'Looks good on the Matlab toolbox side of things! :)\n');
end
fprintf(1,'\n');
input('<<Press any key to continue>>')

% ------------------------------------------------------------------------------
%% 3. Attempt to compile the executables required by the periphery Toolboxes:
% ------------------------------------------------------------------------------
fprintf(1,['\n-2- Compile the binary executables needed for evaluating ' ...
                                                'some operations.\n']);
fprintf(1,['Please make sure that mex is set up with the right compilers for' ...
                                                            ' this system.\n']);
fprintf(1,['Note that errors here are not the end of the world,\nbut mean that ' ...
                        'some operations may fail to execute correctly...\n']);
input('<<Press any key to continue>>')
cd Toolboxes
compile_mex
cd('../');

%-------------------------------------------------------------------------------
%
fprintf(1,'Hope everything compiled ok?!\n\n');

fprintf(1,['All done! Ready when you are to initiate hctsa analysis\nusing a time-series dataset: ' ...
                            'e.g.: TS_Init(''INP_test_ts.mat'')\n']);

% Attempt to add a time series
% SQL_Add('ts','INP_test_ts.txt')
