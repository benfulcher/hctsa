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

	cd OpenTSTOOL/mex-dev
	makemex
	cd ../../
	cd MSmall_utilities
	mex complexitybs.c
	mex nearest.c
	mex shannon.c
	cd ../gpml
    % mex solve_chol.c -largeArrayDims -lmwlapack % !!!! This one won't compile!
	mex sq_dist.c
	cd ../

	cd steps_bumps_toolkit
	mex kvsteps_core.cpp
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
		fprintf(1,'None of the rp_plugins worked. So you will be unable to use this metric\n')
	end

function status = rp_failed()
	plugin_path = fileparts(which('rp_plugin'));

	if ispc
	    plugin_name = 'rp_plugin.exe';
	else
	    plugin_name = 'rp_plugin';
	    if isunix, fileattrib([plugin_path,filesep,plugin_name], '+x'), end
	end

    try
          [status plugin_version] = system([plugin_path,filesep,plugin_name,' -V']);
    catch
          status = 1;
    end
end
end