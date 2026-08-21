function [status, res] = BF_TiseanSystem(tiseanCommand, timeoutSecs)
% BF_TiseanSystem   Run a TISEAN command-line call with a bounded timeout.
%
% All TISEAN-based operations shell out to a compiled binary via system().
% Verified empirically that these binaries do not always error on a bad
% invocation: false_nearest given a missing/unreadable input file blocks
% indefinitely (likely falling back to reading stdin, a common Unix CLI
% pattern, rather than erroring) instead of exiting with a diagnosable
% failure. Left unguarded, that silently stalls a batch computation forever
% with no error, warning, or partial output -- worse than any exception.
%
% This wraps the call with the 'timeout' (or 'gtimeout', GNU coreutils'
% macOS/Homebrew name for the same utility, since it clashes with no
% built-in on that platform) shell command, so a hang fails loudly and
% boundedly instead. If neither is available, degrades to a plain
% unbounded system() call -- identical to every caller's previous
% behavior -- rather than hard-requiring a new dependency.
%
% The default (600s) is deliberately generous, not tight: most of the
% Operations calling this have no length cap on their input (unlike
% NL_RQA/NL_LyapSpec's explicit maxN), so a legitimately long series can
% legitimately take a while, and a timeout that fires on a merely-slow (not
% hung) call would be worse than the bug this fixes -- it would silently
% truncate mid-computation rather than erroring outright (see below).
% Measured directly: NL_c1, the slowest of the 16 TISEAN-calling operations
% tested, took 38s at N=10000 (already longer than any series in this
% toolbox's two standard validation datasets) and 88s at N=30000 -- cost
% flattens with N rather than exploding, so 600s leaves a wide margin
% without being effectively unbounded.
%
% A call that fails this way -- status 124 (this wrapper's own timeout),
% or 126/127 (the shell's reserved "found but not executable"/"command not
% found" codes, i.e. TISEAN never ran at all) -- errors immediately rather
% than returning res to the caller. This is deliberate: a caller expecting
% to parse numeric rows out of res could otherwise misinterpret a partially
% flushed line (captured right before a timeout kill) as real, if
% incomplete, TISEAN output, producing a silently wrong feature value
% instead of a diagnosable failure. Any other exit code means the TISEAN
% binary itself ran to completion (whether cleanly or with its own
% terse refusal, e.g. false_nearest's exit 54 for "too large" -- see
% NL_FNN.m) and is passed through as before for the caller to interpret.
%
%---INPUTS:
% tiseanCommand, the TISEAN command-line string to run (as previously
%                 passed directly to system())
% timeoutSecs [opt], wall-clock time limit in seconds (default: 600)
%
%---OUTPUTS:
% status, the exit status, as returned by system() for the TISEAN binary's
%          own exit (this function already errors on 124/126/127, so a
%          caller only ever sees a status the binary itself chose to exit
%          with).
% res, the captured stdout+stderr text.

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

persistent timeoutCmd
if isempty(timeoutCmd)
    % Cache which (if either) timeout utility is on the path -- checked
    % once per MATLAB session, not on every call:
    if system('command -v timeout >/dev/null 2>&1') == 0
        timeoutCmd = 'timeout';
    elseif system('command -v gtimeout >/dev/null 2>&1') == 0
        timeoutCmd = 'gtimeout';
    else
        timeoutCmd = ''; % neither available -- degrade to unbounded system()
    end
end

if nargin < 2 || isempty(timeoutSecs)
    timeoutSecs = 600;
end

if isempty(timeoutCmd)
    fullCommand = tiseanCommand;
else
    fullCommand = sprintf('%s %g %s', timeoutCmd, timeoutSecs, tiseanCommand);
end

clearLdLibraryPath = isunix && ~ismac;
if clearLdLibraryPath
    % Same fix as BF_RipserSystem.m: MATLAB on Linux prepends its own
    % bundled (older) shared libraries (libgcc_s, libgfortran, etc.) to
    % LD_LIBRARY_PATH, which can shadow the system versions for any
    % subprocess launched via system() -- these TISEAN binaries, freshly
    % compiled by the system C/Fortran compilers, can then fail to load at
    % runtime (running but producing no output) before ever touching the
    % input series. Clearing LD_LIBRARY_PATH just for this subprocess lets
    % the dynamic linker fall back to the normal system search paths. Not
    % needed on macOS, which uses DYLD_LIBRARY_PATH and doesn't hit this.
    %
    % Set via setenv/getenv rather than shell syntax like 'VAR= cmd' -- that
    % inline-assignment form is sh/bash-only and errors under csh/tcsh
    % (which MATLAB's system() uses when $SHELL is csh-family), same as
    % compile_tisean.m already documents for its CC/FC handling.
    oldLdLibraryPath = getenv('LD_LIBRARY_PATH');
    setenv('LD_LIBRARY_PATH','');
end

[status, res] = system(fullCommand);

if clearLdLibraryPath
    setenv('LD_LIBRARY_PATH', oldLdLibraryPath);
end

% See the header comment: any of these three means TISEAN never produced a
% real answer, so error() immediately rather than letting a caller's own
% parsing logic risk mistaking truncated/absent output for a real one.
if any(status == [124, 126, 127])
    switch status
        case 124
            reason = sprintf('timed out after %g seconds', timeoutSecs);
        case 126
            reason = 'found but not executable';
        case 127
            reason = 'not found on the system path -- is TISEAN installed and compiled? (see install.m)';
    end
    error('BF_TiseanSystem:tiseanFailed', ...
        'TISEAN command could not be completed (%s).\nCommand: %s\nOutput: %s', ...
        reason, tiseanCommand, res);
end

end
