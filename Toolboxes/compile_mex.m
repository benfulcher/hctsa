% This helper script compiles the files
% RUN IT FROM Toolboxes!!
% Romesh Abeysuriya, March 2013
% Modified by Ben Fulcher, 2013

tooldir = [pwd,'/'];
fprintf('Now compiling the toolboxes...\n')
if isempty(regexp(tooldir,'Toolboxes'))
    error('This function must be run in the ''Toolboxes'' directory of the HCTSA package...')
end
% fprintf('(I hope %s is the ''Toolboxes'' directory or we have a problem)\n', tooldir)

% Max Little's fastdfa code
fprintf(1,'fastdfa...');
try
    cd([tooldir, 'Max_Little/fastdfa']);
	mex ML_fastdfa_core.c
    fprintf(1,' done.\n');
catch
	error('\nAn error occurred while compiling. Get ''mex ML_fastdfa_core.c'' to work, and then re-run compile.m');
end

% Max Little's Steps Bumps Toolkit
fprintf(1,'Max Little''s ''Steps and bumps'' toolkit...');
cd([tooldir,'Max_Little/steps_bumps_toolkit'])
anyerrors = 0;
try
    mex ML_kvsteps_core.cpp
catch
    fprintf(1,'ERROR: Max Little''s ''Steps and bumps'' toolkit failed to compile correctly\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% Michael Small's code
fprintf(1,'Michael Small''s code...');
cd([tooldir,'Michael_Small'])
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

% Gaussian Process code, gpml
fprintf(1,'Gaussian Process Toolbox, Carl Edward Rasmussen and Hannes Nickisch...');
cd([tooldir,'gpml/util'])
anyerrors = 0;
try
    make
catch
    fprintf(1,'ERROR: Gaussian Process Toolbox failed to compile correctly\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% TSTOOL routines (such a mess)
fprintf(1,'TSTOOL...');
cd([tooldir,'OpenTSTOOL/mex-dev'])
anyerrors = 0;
try
    makemex
catch
    fprintf(1,'ERROR: TSTOOL failed to compile correctly (not a huge surprise)\n');
end
if ~anyerrors, fprintf(1,' done.\n'); end

% TISEAN
% fprintf(1,'Attempting to install TISEAN...!\n');
% cd([tooldir,'Tisean_3.0.1'])
% system('./configure')
% system('make')
% system('make install')


% Return to base directory
cd(tooldir);

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