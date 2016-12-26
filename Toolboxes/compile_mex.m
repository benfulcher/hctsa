% compile_mex   Compiles mex files required for hctsa package
%
% This script must be run in the Toolboxes directory.
%
% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check we're in the correct directory
% ------------------------------------------------------------------------------

currentDir = pwd;

% Path split using platform-dependent separator
weHere = regexp(currentDir,filesep,'split');

if ~strcmp(weHere{end},'Toolboxes')
    error('This code must be run in the ''Toolboxes'' directory of the HCTSA package...')
end

% Sweet. Toolbox path is:
toolDir = [currentDir,filesep];

% ------------------------------------------------------------------------------
% Max Little's fastdfa code
% ------------------------------------------------------------------------------
fprintf(1,'fastdfa...');
try
    cd([toolDir,'Max_Little',filesep,'fastdfa']);
	mex ML_fastdfa_core.c
    fprintf(1,' done.\n');
catch emsg
    fprintf(1,'%s\n\n',emsg.message);
    cd(toolDir);
    errMsg = sprintf(['An error occurred while compiling ML_Fastdfa_core C code.\n' ...
        'It appears that mex is not set up to work on this system (cf. ''doc mex'' and ''mex -setup'').\n' ...
        'Get ''mex ML_fastdfa_core.c'' to work, and then re-run compile_mex.m']);
    error(errMsg)
end

% ------------------------------------------------------------------------------
% Max Little's Steps Bumps Toolkit
% ------------------------------------------------------------------------------
fprintf(1,'Max Little''s ''Steps and bumps'' toolkit...');
cd([toolDir,'Max_Little',filesep,'steps_bumps_toolkit'])
anyErrors = 0;
try
    mex ML_kvsteps_core.cpp
catch
    fprintf(1,'ERROR: Max Little''s ''Steps and bumps'' C++ code failed to compile correctly\n');
end
if ~anyErrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% Max Little's RPDE toolkit
% ------------------------------------------------------------------------------
fprintf(1,'Max Little''s ''RPDE'' code...');
cd([toolDir,'Max_Little',filesep,'rpde'])
anyErrors = 0;
try
    mex ML_close_ret.c
catch
    fprintf(1,'ERROR: Max Little''s ''RPDE'' C code failed to compile correctly\n');
end
if ~anyErrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% Michael Small's code
% ------------------------------------------------------------------------------
fprintf(1,'Michael Small''s code...');
cd([toolDir,'Michael_Small'])
anyErrors = 0;
try
    mex MS_complexitybs.c % compile Michael Small's complexitybs C code
catch
    fprintf(1,'ERROR: Michael Small''s ''complexitybs'' C code failed to compile correctly\n');
end
try
    mex MS_nearest.c      % compile Michael Small's nearest C code
catch
    fprintf(1,'ERROR: Michael Small''s ''nearest'' C code failed to compile correctly\n');
end
try
    mex MS_shannon.c      % compile Michael Small's shannon C code
catch
    fprintf(1,'ERROR: Michael Small''s ''shannon'' C code failed to compile correctly\n');
end
if ~anyErrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% Gaussian Process code, gpml
% ------------------------------------------------------------------------------
fprintf(1,'Gaussian Process Toolbox, Carl Edward Rasmussen and Hannes Nickisch...');
cd([toolDir,'gpml',filesep,'util'])
anyErrors = 0;
try
    make
catch
    fprintf(1,'ERROR: Gaussian Process Toolbox failed to compile correctly\n');
end
if ~anyErrors, fprintf(1,' done.\n'); end

%-------------------------------------------------------------------------------
% Physionet sample entropy code (turned to mex)
%-------------------------------------------------------------------------------
fprintf(1,'Sample entropy...');
cd([toolDir,'Physionet'])
anyErrors = 0;
try
    mex sampen_mex.c
catch
    fprintf(1,'ERROR: Physionet implementation of sample entropy failed to compile.\n');
end
if ~anyErrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% TSTOOL routines (such a mess)
% ------------------------------------------------------------------------------
fprintf(1,'TSTOOL...');
cd([toolDir,'OpenTSTOOL',filesep,'mex-dev'])
anyErrors = 0;
try
    makemex
    settspath(fullfile(toolDir,'OpenTSTOOL'));
catch
    fprintf(1,'ERROR: The TSTOOL package failed to compile correctly\n');
end
if ~anyErrors, fprintf(1,'Done!\n'); end

% ------------------------------------------------------------------------------
% TISEAN
% ------------------------------------------------------------------------------
fprintf(1,['NB: To use TISEAN routines, you need to compile and install them on your system\n' ...
    'In the commandline, navigate to the Toolboxes%sTisean_3.0.1 directory\n  ''./configure'',' ...
    ' ''make'' and ''make install'' commands\n(cf. Documentation for instructions,' ...
    ' including implementations for Windows)...\n'],filesep);

% Return to base directory
cd(toolDir);
