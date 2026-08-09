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
fprintf(1,'TISEAN nonlinear time-series analysis routines... (not mex, compilation, be patient...)');
tiseanDir = fullfile(toolDir,'Tisean_3.0.1');
cd(tiseanDir);
ok = true;
if ispc
    ok = false;
    fprintf(1,['\nERROR: automatic compilation of TISEAN is not supported on Windows.\n' ...
        'See Toolboxes%sTisean_3.0.1%sindex.html for manual build instructions.\n'],filesep,filesep);
elseif ~compilerOnPath('gcc') || ~compilerOnPath('make')
    % Check for a C compiler and make up front, rather than letting TISEAN's
    % ./configure fail and dump a wall of cryptic autoconf output.
    % NB: don't check via 'which'/'command -v' + shell redirection -- MATLAB's
    % system() runs the command through $SHELL, and on clusters that's often
    % csh/tcsh, which understands neither sh's '2>&1' redirect syntax ("Ambiguous
    % output redirect") nor the 'command -v' builtin (sh/bash-only). Both give a
    % false "no compiler found" even when gcc/make are genuinely on PATH.
    % Instead, actually try to run the compiler and check the exit status;
    % suppress its output via system()'s own second return value (captured at
    % the MATLAB level, not via shell redirection), which works under any shell.
    ok = false;
    if ismac
        fprintf(1,['\nERROR: no C compiler (gcc/clang) or ''make'' found on the system PATH.\n' ...
            'Install the Xcode Command Line Tools with: xcode-select --install\n' ...
            'then re-run compile_mex.m.\n']);
    else
        fprintf(1,['\nERROR: no C compiler (gcc) or ''make'' found on the system PATH.\n' ...
            'On Debian/Ubuntu: sudo apt-get install build-essential\n' ...
            'then re-run compile_mex.m.\n']);
    end
    % If gcc/make are on PATH in an interactive terminal but not here, MATLAB's
    % subprocess environment differs from your shell's -- common on clusters when
    % MATLAB runs on a different node (e.g., a batch job landed on a compute node
    % lacking build tools) or a job scheduler didn't export your full environment.
    % Print what MATLAB itself sees so that mismatch is diagnosable without
    % needing to reproduce it interactively:
    [~,hostnameOut] = system('hostname');
    [~,pathOut] = system('echo $PATH');
    fprintf(1,['\nDiagnostics (compare against your interactive terminal on the same node):\n' ...
        '  hostname: %s' ...
        '  PATH:     %s\n'],hostnameOut,pathOut);
else
    % Build against a prefix outside the repo, in case the repo lives under a
    % path containing spaces (e.g., a Dropbox-synced folder): TISEAN's 1990s-era
    % configure/Makefiles/install-sh don't quote paths internally, so an
    % install prefix with spaces breaks 'make install'. Compiling itself is
    % unaffected by spaces in the source path, so we build in place and only
    % route the *install* step through a space-free temporary prefix, then
    % copy the resulting binaries into place with MATLAB's own copyfile.
    buildPrefix = tempname;
    mkdir(fullfile(buildPrefix,'bin'));
    % Modern C compilers (e.g., Xcode Clang) reject this package's K&R-style
    % implicit-int/implicit-function-declaration code by default; -std=gnu89
    % restores the old C89 semantics it was written against.
    % Set CC via MATLAB's own setenv rather than shell syntax like 'CC="..." cmd'
    % -- that inline-assignment form is sh/bash-only and errors under csh/tcsh
    % (which MATLAB's system() uses when $SHELL is csh-family), since csh has no
    % such syntax and instead tries to run "CC=..." as a literal command name.
    % setenv() sets the process environment directly, so configure (itself a
    % #!/bin/sh script) inherits CC correctly regardless of the caller's shell.
    oldCC = getenv('CC');
    tiseanCC = 'gcc -std=gnu89 -Wno-implicit-int -Wno-implicit-function-declaration';
    setenv('CC',tiseanCC);
    [statusConfigure,outConfigure] = system(sprintf('./configure --prefix=%s',buildPrefix));
    setenv('CC',oldCC);
    if statusConfigure~=0
        ok = false;
        fprintf(1,'\nERROR: TISEAN ./configure failed:\n%s\n',outConfigure);
    else
        [statusMake,outMake] = system('make');
        if statusMake~=0
            ok = false;
            fprintf(1,'\nERROR: TISEAN make failed:\n%s\n',outMake);
        else
            % NB: 'make install' runs its own "missing" self-check that
            % reports false positives (it probes each binary by bare name
            % before bin/ is on PATH), so verify success directly against
            % buildPrefix/bin instead of trusting its exit status or missing.log.
            system('make install');
            keyBinaries = {'c1','d2','c2d','c2g','c2t','nstat_z','false_nearest'};
            builtOk = all(cellfun(@(f)isfile(fullfile(buildPrefix,'bin',f)),keyBinaries));
            if builtOk
                if ~isfolder('bin')
                    mkdir('bin');
                end
                copyfile(fullfile(buildPrefix,'bin'),'bin');
                system('chmod +x bin/*');
                fprintf(1,' done (binaries installed to %s).\n',fullfile(tiseanDir,'bin'));
            else
                ok = false;
                fprintf(1,'\nERROR: TISEAN build finished but required binaries are missing.\n');
            end
        end
    end
    rmdir(buildPrefix,'s');
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

function tf = compilerOnPath(cmdName)
    % Checks whether cmdName is actually invocable, without relying on
    % 'which'/'command -v' or shell redirection syntax -- both differ between
    % sh and csh/tcsh, and MATLAB's system() uses $SHELL. Just run the tool
    % with --version and check the exit status; the two-output form of
    % system() captures (and suppresses) output at the MATLAB level, so no
    % shell-specific redirection is needed.
    [status,~] = system([cmdName,' --version']);
    tf = (status==0);
end
