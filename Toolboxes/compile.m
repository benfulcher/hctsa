% This helper script compiles the files
% RUN IT FROM Toolboxes!!
% Romesh Abeysuriya, March 2013
% function compile

tooldir = [pwd,'/'];
fprintf('Now compiling the toolboxes...\n')
fprintf('(I hope %s is the Toolbox directory or we have a problem)\n', tooldir)

% Max Little's fastdfa code
fprintf(1,'fastdfa...');
try
    cd([tooldir, 'Max_Little/fastdfa']);
	mex fastdfa_core.c
    fprintf(1,' done.\n');
catch
	error('\nAn error occurred while compiling. Get ''mex fastdfa_core.c'' to work, and then re-run compile.m');
end

% MICHAEL SMALL
fprintf(1,'Michael Small''s code...');
cd([tooldir,'Michael_Small'])
mex MS_complexitybs.c % compile Michael Small's complexitybs C code
mex MS_nearest.c % compile Michael Small's nearest C code
mex MS_shannon.c % compile Michael Small's shannon C code
fprintf(1,' done.\n');

% Gaussian Process code, gpml
fprintf(1,'Gaussian Process Toolbox, Carl Edward Rasmussen & Hannes Nickisch...');
cd([tooldir,'gpml/util'])
make
fprintf(1,' done.\n');

% Max Little's Steps Bumps Toolkit
fprintf(1,'Steps and bumps toolkit, Max Little...');
cd([tooldir,'steps_bumps_toolkit'])
mex ML_kvsteps_core.cpp
fprintf(1,' done.\n');

% TSTOOL
fprintf(1,'TSTOOL...');
cd([tooldir,'OpenTSTOOL/mex-dev'])
makemex
fprintf(1,' done.\n');

    
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