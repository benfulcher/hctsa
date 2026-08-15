function ok = compile_tisean()
% compile_tisean   Compiles and installs the TISEAN nonlinear time-series toolkit.
%
% This script must be run in the Toolboxes directory (same requirement as
% compile_mex.m, which calls this as part of its own build sequence).
%
% Split out from compile_mex.m into its own function so that, if TISEAN is
% the only component that failed (e.g. a required compiler was missing),
% you can fix that and just re-run this function on its own -- rather than
% re-running (and re-compiling) every other component compile_mex.m
% handles.
%
% Only builds the specific TISEAN binaries hctsa actually calls (grepped
% from Operations/*.m), not the full ~70-tool TISEAN suite -- see the
% comment above the build step for why.
%
% ---OUTPUTS:
% ok, true if TISEAN compiled and all binaries hctsa depends on (c1, c2d,
%     c2g, c2t, d2, nstat_z, false_nearest, boxcount, lyap_r, lyap_spec,
%     poincare, nrlazy) were installed.

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
weHere = regexp(currentDir,filesep,'split');
if ~strcmp(weHere{end},'Toolboxes')
    error('This code must be run in the ''Toolboxes'' directory of the HCTSA package...')
end
toolDir = currentDir;

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
            'then re-run compile_tisean.m.\n']);
    else
        fprintf(1,['\nERROR: no C compiler (gcc) or ''make'' found on the system PATH.\n' ...
            'On Debian/Ubuntu: sudo apt-get install build-essential\n' ...
            'then re-run compile_tisean.m.\n']);
    end
    printPathDiagnostics();
elseif isempty(findFortranCompiler())
    % Also required: 4 of the 12 binaries hctsa depends on (c1, c2d, c2g, c2t)
    % are Fortran, not C (TISEAN's source_f/, not source_c/), so a C-only
    % toolchain isn't sufficient for a complete build.
    ok = false;
    if ismac
        fprintf(1,['\nERROR: no Fortran compiler (gfortran/g77) found on the system PATH ' ...
            '(also checked the standard Homebrew install locations directly).\n' ...
            'TISEAN needs this in addition to a C compiler: 4 of the 12 binaries hctsa\n' ...
            'depends on (c1, c2d, c2g, c2t) are Fortran, not C.\n' ...
            'Install gfortran (e.g. via Homebrew: brew install gcc), then re-run compile_tisean.m.\n' ...
            'If it''s already installed via Homebrew but still not found: MATLAB launched from\n' ...
            'Finder/Dock does not inherit your shell profile''s PATH, which is where Homebrew adds\n' ...
            'itself -- this check already looks in /opt/homebrew/bin and /usr/local/bin directly to\n' ...
            'cover that, so if it''s still not found, your gfortran may be somewhere else.\n']);
    else
        fprintf(1,['\nERROR: no Fortran compiler (gfortran/g77) found on the system PATH.\n' ...
            'TISEAN needs this in addition to a C compiler: 4 of the 12 binaries hctsa\n' ...
            'depends on (c1, c2d, c2g, c2t) are Fortran, not C.\n' ...
            'On Debian/Ubuntu: sudo apt-get install gfortran\n' ...
            'then re-run compile_tisean.m.\n']);
    end
    printPathDiagnostics();
else
    % Only build the specific TISEAN binaries hctsa actually calls (grepped
    % from Operations/*.m: 4 Fortran, 8 C), rather than all ~70 tools TISEAN
    % ships. This isn't just faster -- most of TISEAN's other tools are
    % 1990s-era code that produces (mostly harmless, legacy-syntax) compiler
    % warnings under modern compilers, and one (cluster.f) outright fails to
    % compile on modern gfortran (a bad Fortran format string that older
    % gfortran tolerated as an extension but current versions reject). None
    % of that matters if those files are simply never compiled.
    fortranBins = {'c1','c2d','c2g','c2t'};
    cBins = {'d2','nstat_z','false_nearest','boxcount','lyap_r','lyap_spec','poincare','nrlazy'};

    % Modern C compilers (e.g., Xcode Clang) reject this package's K&R-style
    % implicit-int/implicit-function-declaration code by default; -std=gnu89
    % restores the old C89 semantics it was written against. -std=legacy is
    % the Fortran analogue: without it, modern gfortran emits a warning for
    % every old-style shared/labeled DO-loop termination in this 1990s F77
    % code (a large but functionally harmless volume, given -std=legacy
    % still fully supports these constructs as a GNU extension -- it isn't
    % suppressing anything at risk of being silently wrong, just silencing
    % pedantry about outdated-but-supported syntax).
    % Set CC/FC via MATLAB's own setenv rather than shell syntax like
    % 'CC="..." cmd' -- that inline-assignment form is sh/bash-only and
    % errors under csh/tcsh (which MATLAB's system() uses when $SHELL is
    % csh-family), since csh has no such syntax and instead tries to run
    % "CC=..." as a literal command name. setenv() sets the process
    % environment directly, so configure (itself a #!/bin/sh script)
    % inherits CC/FC correctly regardless of the caller's shell.
    oldCC = getenv('CC');
    oldFC = getenv('FC');
    tiseanCC = 'gcc -std=gnu89 -Wno-implicit-int -Wno-implicit-function-declaration';
    tiseanFC = [findFortranCompiler(),' -std=legacy']; % may be an absolute path (see findFortranCompiler)
    setenv('CC',tiseanCC);
    setenv('FC',tiseanFC);
    % Clear any cached configure results from a previous run (e.g. one that
    % ran before a missing compiler was installed) so CC/FC get freshly
    % re-detected rather than reusing a stale cache:
    if isfile('config.cache')
        delete('config.cache');
    end
    % --prefix is required by configure but never actually used below (we
    % copy binaries directly from source_c//source_f/ ourselves, bypassing
    % 'make install' entirely -- see below), so a nonexistent placeholder is
    % fine:
    [statusConfigure,outConfigure] = system(sprintf('./configure --prefix=%s',tempname));
    setenv('CC',oldCC);
    setenv('FC',oldFC);
    if statusConfigure~=0
        ok = false;
        fprintf(1,'\nERROR: TISEAN ./configure failed:\n%s\n',outConfigure);
    else
        [~,outMakeF] = system(['make -C source_f ',strjoin(fortranBins,' ')]);
        [~,outMakeC] = system(['make -C source_c ',strjoin(cBins,' ')]);
        % NB: TISEAN's own per-binary compile recipes are '-'-prefixed (make
        % ignores a failed individual compile and continues), so neither
        % make's exit status nor missing.log-style self-checks are reliable
        % -- verify success directly by checking each binary actually exists
        % where it was just built:
        keyBinaries = [fortranBins,cBins];
        builtPath = [cellfun(@(f)fullfile('source_f',f),fortranBins,'UniformOutput',false), ...
            cellfun(@(f)fullfile('source_c',f),cBins,'UniformOutput',false)];
        gotBinary = cellfun(@isfile,builtPath);
        if all(gotBinary)
            if ~isfolder('bin')
                mkdir('bin');
            end
            for i = 1:numel(builtPath)
                copyfile(builtPath{i},fullfile('bin',keyBinaries{i}));
            end
            system('chmod +x bin/*');
            fprintf(1,' done\n(binaries installed to %s).\n',fullfile(tiseanDir,'bin'));
        else
            ok = false;
            fprintf(1,'\nERROR: TISEAN build finished but required binaries are missing: %s\n%s%s', ...
                strjoin(keyBinaries(~gotBinary),', '),outMakeF,outMakeC);
        end
    end
end
cd(toolDir);

end

% ------------------------------------------------------------------------------
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

% ------------------------------------------------------------------------------
function fc = findFortranCompiler()
    % Looks for gfortran, then g77: first via PATH (as gcc/make are checked),
    % then -- on macOS only -- directly in the two standard Homebrew install
    % prefixes (/opt/homebrew for Apple Silicon, /usr/local for Intel).
    %
    % The Homebrew-prefix fallback covers a common, well-documented macOS
    % gotcha, not a machine-specific assumption: MATLAB launched from
    % Finder/Dock (rather than a terminal) runs outside a login shell, so it
    % never sources the PATH additions Homebrew's installer appends to
    % ~/.zprofile -- a genuinely-installed Homebrew gfortran can therefore be
    % invisible to MATLAB's system() even though `gfortran --version` works
    % fine in an interactive terminal. This isn't needed for gcc/make, which
    % come from Xcode Command Line Tools and install to /usr/bin -- part of
    % the OS-level default PATH regardless of how the process was launched.
    fc = '';
    candidates = {'gfortran','g77'};
    for i = 1:numel(candidates)
        if compilerOnPath(candidates{i})
            fc = candidates{i};
            return
        end
    end
    if ismac
        homebrewPrefixes = {'/opt/homebrew/bin','/usr/local/bin'};
        for i = 1:numel(homebrewPrefixes)
            for j = 1:numel(candidates)
                candidatePath = fullfile(homebrewPrefixes{i},candidates{j});
                if isfile(candidatePath)
                    [status,~] = system(['"',candidatePath,'" --version']);
                    if status==0
                        fc = candidatePath;
                        return
                    end
                end
            end
        end
    end
end

% ------------------------------------------------------------------------------
function printPathDiagnostics()
    % If a compiler is on PATH in an interactive terminal but not here,
    % MATLAB's subprocess environment differs from your shell's -- common on
    % clusters when MATLAB runs on a different node (e.g., a batch job
    % landed on a compute node lacking build tools) or a job scheduler
    % didn't export your full environment. Print what MATLAB itself sees so
    % that mismatch is diagnosable without needing to reproduce it
    % interactively:
    [~,hostnameOut] = system('hostname');
    [~,pathOut] = system('echo $PATH');
    fprintf(1,['\nDiagnostics (compare against your interactive terminal on the same node):\n' ...
        '  hostname: %s' ...
        '  PATH:     %s\n'],hostnameOut,pathOut);
end
