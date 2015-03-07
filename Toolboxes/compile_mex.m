% This script compiles the mex files required for all operations
% implemented in the HCTSA package.
% It must be run in the Toolboxes directory.
% 
%---HISTORY:
% Tweaks by Dror Cohen, 2014-04-08
% Modified by Ben Fulcher, 2013
% Romesh Abeysuriya, March 2013
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% Check we're in the correct folder
% ------------------------------------------------------------------------------

currentDir = pwd;

% Path split using platform-dependent separator
if isunix
    weHere = regexp(currentDir,'/','split');
else
    weHere = regexp(currentDir,'\','split');
end

if ~strcmp(weHere{end},'Toolboxes')
    error('This code must be run in the ''Toolboxes'' directory of the HCTSA package...')
end

% Sweet. Toolbox path is:
toolDir = [currentDir '/'];

% ------------------------------------------------------------------------------
% Max Little's fastdfa code
% ------------------------------------------------------------------------------
fprintf(1,'fastdfa...');
try
    cd([toolDir, 'Max_Little/fastdfa']);
	mex ML_fastdfa_core.c
    fprintf(1,' done.\n');
catch emsg
    fprintf(1,'%s\n\n',emsg.message);
    error(['An error occurred while compiling.\n' ...
        'It appears that mex is not set up to work on this system (cf. ''doc mex'' and ''mex -setup'').\n' ...
        'Get ''mex ML_fastdfa_core.c'' to work, and then re-run compile_mex.m']);
end

% ------------------------------------------------------------------------------
% Max Little's Steps Bumps Toolkit
% ------------------------------------------------------------------------------
fprintf(1,'Max Little''s ''Steps and bumps'' toolkit...');
cd([toolDir,'Max_Little/steps_bumps_toolkit'])
anyerrors = 0;
try
    mex ML_kvsteps_core.cpp
catch
    fprintf(1,'ERROR: Max Little''s ''Steps and bumps'' toolkit failed to compile correctly\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% Michael Small's code
% ------------------------------------------------------------------------------
fprintf(1,'Michael Small''s code...');
cd([toolDir,'Michael_Small'])
anyerrors = 0;
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
if ~anyerrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% Gaussian Process code, gpml
% ------------------------------------------------------------------------------
fprintf(1,'Gaussian Process Toolbox, Carl Edward Rasmussen and Hannes Nickisch...');
cd([toolDir,'gpml/util'])
anyerrors = 0;
try
    make
catch
    fprintf(1,'ERROR: Gaussian Process Toolbox failed to compile correctly\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% TSTOOL routines (such a mess)
% ------------------------------------------------------------------------------
fprintf(1,'TSTOOL...');
cd([toolDir,'OpenTSTOOL/mex-dev'])
anyerrors = 0;
try
    makemex
catch
    fprintf(1,'ERROR: TSTOOL failed to compile correctly (not a huge surprise)\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% ------------------------------------------------------------------------------
% TISEAN
% ------------------------------------------------------------------------------
fprintf(1,['NB: To use TISEAN routines, you need to install them on your system using\n  ''./configure''' ...
            ' ''make clean'' and ''make'' and ''make install'' commands (cf. Documentation)...\n']);

% Return to base directory
cd(toolDir);

%     fprintf('------\nNow installing the CRPTool\n');
%     fprintf('At the prompts\n-Remove existing toolbox\n-Create folder if asked\n-DO NOT add to the path\n-DO NOT delete the installation file\n-----\n');
% 
%     plugininstall_x86_64('./CRPTool');
%     if rp_failed()
%         fprintf(1,'rp_plugin failed, trying a different version...\n')
%         plugininstall_i686('./CRPTool');
%     else 
%         return
%     end
%     if rp_failed()
%         fprintf(1,'rp_plugin failed, trying a different version...\n')
%         plugininstall_ia64('./CRPTool');
%     else 
%         return
%     end
%     if rp_failed()
%         fprintf(1,'None of the rp_plugins worked. So you will be unable to use this operation\n')
%     end
% 
% function status = rp_failed()
%     plugin_path = fileparts(which('rp_plugin'));
% 
%     if ispc
%         plugin_name = 'rp_plugin.exe';
%     else
%         plugin_name = 'rp_plugin';
%         if isunix
%             fileattrib([plugin_path, filesep, plugin_name], '+x')
%         end
%     end
% 
%     try
%           [status, plugin_version] = system([plugin_path,filesep,plugin_name,' -V']);
%     catch
%           status = 1;
%     end
% end
% end