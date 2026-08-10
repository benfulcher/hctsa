% compile_mex   Compiles mex files required for hctsa package
%
% This script must be run in the Toolboxes directory.
%
% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
toolDir = currentDir;

% Track per-component success so a final summary can be printed, and so
% one early failure (e.g. mex not yet configured) doesn't abort the whole
% script and silently skip every toolbox after it:
results = struct('name',{},'ok',{});

% ------------------------------------------------------------------------------
% Max Little's fastdfa code
% ------------------------------------------------------------------------------
fprintf(1,'fastdfa...');
ok = true;
try
    cd(fullfile(toolDir,'Max_Little','fastdfa'));
	mex ML_fastdfa_core.c
    fprintf(1,' done.\n');
catch emsg
    ok = false;
    fprintf(1,'%s\n',emsg.message);
    fprintf(1,['ERROR: ML_fastdfa_core C code failed to compile. It appears that mex is not ' ...
        'set up to work on this system (cf. ''doc mex'' and ''mex -setup'').\n']);
end
results(end+1) = struct('name','Max Little''s fastdfa','ok',ok);

% ------------------------------------------------------------------------------
% Max Little's Steps Bumps Toolkit
% ------------------------------------------------------------------------------
fprintf(1,'Max Little''s ''Steps and bumps'' toolkit...');
cd(fullfile(toolDir,'Max_Little','steps_bumps_toolkit'))
ok = true;
try
    mex ML_kvsteps_core.cpp
    fprintf(1,' done.\n');
catch
    ok = false;
    fprintf(1,'ERROR: Max Little''s ''Steps and bumps'' C++ code failed to compile correctly.\n');
end
results(end+1) = struct('name','Max Little''s steps_bumps toolkit','ok',ok);

% ------------------------------------------------------------------------------
% Max Little's RPDE toolkit
% ------------------------------------------------------------------------------
fprintf(1,'Max Little''s ''RPDE'' code...');
cd(fullfile(toolDir,'Max_Little','rpde'))
ok = true;
try
    mex ML_close_ret.c
    fprintf(1,' done.\n');
catch
    ok = false;
    fprintf(1,'ERROR: Max Little''s ''RPDE'' C code failed to compile correctly\n');
end
results(end+1) = struct('name','Max Little''s RPDE toolkit','ok',ok);

% ------------------------------------------------------------------------------
% Michael Small's code
% ------------------------------------------------------------------------------
fprintf(1,'Michael Small''s code...');
cd(fullfile(toolDir,'Michael_Small'))
ok = true;
try
    mex MS_complexitybs.c % compile Michael Small's complexitybs C code
catch
    ok = false;
    fprintf(1,'ERROR: Michael Small''s ''complexitybs'' C code failed to compile correctly\n');
end
try
    mex MS_nearest.c      % compile Michael Small's nearest C code
catch
    ok = false;
    fprintf(1,'ERROR: Michael Small''s ''nearest'' C code failed to compile correctly\n');
end
try
    mex MS_shannon.c      % compile Michael Small's shannon C code
catch
    ok = false;
    fprintf(1,'ERROR: Michael Small''s ''shannon'' C code failed to compile correctly\n');
end
if ok
    fprintf(1,' done.\n');
end
results(end+1) = struct('name','Michael Small''s code','ok',ok);

% ------------------------------------------------------------------------------
% Gaussian Process code, gpml
% ------------------------------------------------------------------------------
fprintf(1,'Gaussian Process Toolbox, Carl Edward Rasmussen and Hannes Nickisch...');
cd(fullfile(toolDir,'gpml','util'))
ok = true;
try
    make
    fprintf(1,' done.\n');
catch
    ok = false;
    fprintf(1,'ERROR: Gaussian Process Toolbox failed to compile correctly.\n');
end
results(end+1) = struct('name','Gaussian Process Toolbox (gpml)','ok',ok);

%-------------------------------------------------------------------------------
% Physionet sample entropy code (turned to mex)
%-------------------------------------------------------------------------------
fprintf(1,'Sample entropy...');
cd(fullfile(toolDir,'Physionet'))
ok = true;
try
    mex sampen_mex.c
    fprintf(1,' done.\n');
catch
    ok = false;
    fprintf(1,'ERROR: Physionet implementation of sample entropy failed to compile.\n');
end
results(end+1) = struct('name','Physionet sample entropy','ok',ok);

% ------------------------------------------------------------------------------
% TISEAN
% ------------------------------------------------------------------------------
% (its own function, compile_tisean.m, so it can be re-run on its own if it's
% the only component that fails -- e.g. after installing a missing compiler --
% without re-running every other component here too)
cd(toolDir);
ok = compile_tisean();
if ~ok
    fprintf(1,['If TISEAN was the only component that failed above, fix the issue and just ' ...
        're-run compile_tisean (rather than the whole of compile_mex, which would ' ...
        're-compile everything else too).\n']);
end
results(end+1) = struct('name','TISEAN','ok',ok);

%-------------------------------------------------------------------------------
% CATCH22
%-------------------------------------------------------------------------------
fprintf(1,'catch22...');
cd(fullfile(toolDir,'catch22','wrap_Matlab'))
ok = true;
try
    mexAll
    fprintf(1,' done.\n');
catch emsg
    ok = false;
    fprintf(1,'ERROR: catch22 failed to compile correctly.\n%s\n',emsg.message);
end
results(end+1) = struct('name','catch22','ok',ok);

% Return to base directory
cd(toolDir);

%-------------------------------------------------------------------------------
% Summary
%-------------------------------------------------------------------------------
fprintf(1,'\n---Compilation summary---\n');
for i = 1:length(results)
    if results(i).ok
        fprintf(1,'  [ OK ] %s\n',results(i).name);
    else
        fprintf(1,'  [FAIL] %s\n',results(i).name);
    end
end
numOk = sum([results.ok]);
fprintf(1,'%u/%u components compiled successfully.\n',numOk,length(results));
if numOk < length(results)
    fprintf(1,['Some components failed to compile -- hctsa will still run, but operations that ' ...
        'depend on a failed component will not work. See the ERROR messages above for details.\n']);
end
