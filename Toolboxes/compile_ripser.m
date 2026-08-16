function ok = compile_ripser()
% compile_ripser   Compiles and installs the ripser persistent-homology tool.
%
% This script must be run in the Toolboxes directory (same requirement as
% compile_mex.m, which calls this as part of its own build sequence).
%
% Split out from compile_mex.m into its own function so that, if ripser is
% the only component that failed (e.g. a required compiler was missing), you
% can fix that and just re-run this function on its own -- rather than
% re-running (and re-compiling) every other component compile_mex.m handles.
%
% Unlike TISEAN (see compile_tisean.m), ripser is a single C++ translation
% unit with no ./configure/make/Fortran step -- one compiler invocation.
%
% ---OUTPUTS:
% ok, true if ripser compiled and the 'ripser' binary was installed.

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

fprintf(1,'ripser persistent-homology tool...');
ripserDir = fullfile(toolDir,'ripser');
cd(ripserDir);
ok = true;
if ispc
    ok = false;
    fprintf(1,['\nERROR: automatic compilation of ripser is not supported on Windows.\n' ...
        'See Toolboxes%sripser%sREADME.md for manual build instructions.\n'],filesep,filesep);
else
    cxx = findCXXCompiler();
    if isempty(cxx)
        ok = false;
        if ismac
            fprintf(1,['\nERROR: no C++ compiler (g++/clang++) found on the system PATH.\n' ...
                'Install the Xcode Command Line Tools with: xcode-select --install\n' ...
                'then re-run compile_ripser.m.\n']);
        else
            fprintf(1,['\nERROR: no C++ compiler (g++/clang++) found on the system PATH.\n' ...
                'On Debian/Ubuntu: sudo apt-get install build-essential\n' ...
                'then re-run compile_ripser.m.\n']);
        end
    else
        % Single translation unit -- the exact invocation ripser's own
        % Makefile uses for its default 'ripser' target (no coefficients, no
        % debug symbols):
        binPath = fullfile('bin','ripser');
        if ~isfolder('bin')
            mkdir('bin');
        end
        [statusBuild,outBuild] = system(sprintf('%s -std=c++11 -Wall ripser.cpp -o %s -O3 -D NDEBUG', ...
            cxx,binPath));
        % As with TISEAN, don't trust the compiler's exit status alone --
        % verify the binary actually landed where expected:
        if statusBuild==0 && isfile(binPath)
            system(['chmod +x ',binPath]);
            fprintf(1,' done\n(binary installed to %s).\n',fullfile(ripserDir,'bin'));
        else
            ok = false;
            fprintf(1,'\nERROR: ripser build failed:\n%s\n',outBuild);
        end
    end
end
cd(toolDir);

end

% ------------------------------------------------------------------------------
function tf = compilerOnPath(cmdName)
    % See the identical helper in compile_tisean.m for why this avoids
    % 'which'/'command -v' and shell redirection (csh/tcsh compatibility).
    [status,~] = system([cmdName,' --version']);
    tf = (status==0);
end

% ------------------------------------------------------------------------------
function cxx = findCXXCompiler()
    % Looks for g++, then clang++, on PATH -- and, on macOS only, directly in
    % the standard Homebrew install prefixes as a fallback (see the identical
    % rationale in compile_tisean.m's findFortranCompiler: MATLAB launched
    % from Finder/Dock doesn't inherit a login shell's Homebrew PATH
    % additions). Xcode's own clang++ lives in /usr/bin and needs no such
    % fallback, but a Homebrew-installed g++ might not be found otherwise.
    cxx = '';
    candidates = {'c++','g++','clang++'};
    for i = 1:numel(candidates)
        if compilerOnPath(candidates{i})
            cxx = candidates{i};
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
                        cxx = ['"',candidatePath,'"'];
                        return
                    end
                end
            end
        end
    end
end
