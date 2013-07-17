% This helper script compiles the files
% RUN IT FROM code_sub/Toolboxes!!
% Romesh Abeysuriya, March 2013
function compile
	fprintf('Now compiling the toolboxes\n')
	try
        cd fastdfa
		mex fastdfa_core.c
        cd ../
	catch
		error('An error occurred while compiling. Get ''mex fastdfa_core.c'' to work, and then re-run compile.m');
	end

    % TSTOOL
	cd OpenTSTOOL/mex-dev
	makemex
	cd ../../

    % MICHAEL SMALL
	cd Michael_Small
	mex MS_complexitybs.c % compile Michael Small's complexitybs C code
	mex MS_nearest.c % compile Michael Small's nearest C code
	mex MS_shannon.c % compile Michael Small's shannon C code

    % Gaussian Process code, gpml
	cd ../gpml
	mex sq_dist.c
    % mex solve_chol.c -largeArrayDims -lmwlapack % !!!! This one won't compile!
	cd ../

    % Max Little's Steps Bumps Toolkit
	cd steps_bumps_toolkit
	mex ML_kvsteps_core.cpp
	cd ../

	fprintf('------\nNow installing the CRPTool\n');
	fprintf('At the prompts\n-Remove existing toolbox\n-Create folder if asked\n-DO NOT add to the path\n-DO NOT delete the installation file\n-----\n');

	plugininstall_x86_64('./CRPTool');
	if rp_failed()
		fprintf(1,'rp_plugin failed, trying a different version...\n')
		plugininstall_i686('./CRPTool');
	else 
		return
	end
	if rp_failed()
		fprintf(1,'rp_plugin failed, trying a different version...\n')
		plugininstall_ia64('./CRPTool');
	else 
		return
	end
	if rp_failed()
		fprintf(1,'None of the rp_plugins worked. So you will be unable to use this operation\n')
	end

function status = rp_failed()
	plugin_path = fileparts(which('rp_plugin'));

	if ispc
	    plugin_name = 'rp_plugin.exe';
	else
	    plugin_name = 'rp_plugin';
	    if isunix
            fileattrib([plugin_path,filesep,plugin_name], '+x')
        end
	end

    try
          [status, plugin_version] = system([plugin_path,filesep,plugin_name,' -V']);
    catch
          status = 1;
    end
end
end