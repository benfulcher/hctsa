function [status, res] = BF_TiseanSystem(tiseanCommand, timeoutSecs)
% BF_TiseanSystem   Run a TISEAN command-line call with a bounded timeout.
%
% All TISEAN-based operations shell out to a compiled binary via system().
% These binaries do not always error on a bad invocation: e.g., false_nearest
% given a missing/unreadable input file can block indefinitely instead of
% exiting with a diagnosable failure. Left unguarded, that silently stalls a
% batch computation forever with no error, warning, or partial output --
% worse than any exception.
%
% This wraps the call with the 'timeout' (or 'gtimeout', GNU coreutils'
% macOS/Homebrew name for the same utility, since it clashes with no
% built-in on that platform) shell command, so a hang fails loudly and
% boundedly instead. If neither is available, degrades to a plain
% unbounded system() call (with a one-time warning) rather than
% hard-requiring a new dependency; a post-hoc cap on the captured output
% size then still catches a binary that runs but streams console output
% without terminating.
%
% The default (600s) is deliberately generous, not tight: most of the
% Operations calling this have no length cap on their input (unlike
% NL_RQA/NL_LyapSpec's explicit maxN), so a legitimately long series can
% legitimately take a while, and a timeout that fires on a merely-slow (not
% hung) call would needlessly abort a valid computation (see below). It is a
% cap on computation time for hung calls, not a tight budget.
%
% A call that fails this way -- status 124 (this wrapper's own timeout),
% 126/127 (the shell's reserved "found but not executable"/"command not
% found" codes, i.e. TISEAN never ran at all), or an over-limit output
% capture (a runaway binary, see below) -- errors immediately rather
% than returning res to the caller. This is deliberate: a caller expecting
% to parse numeric rows out of res could otherwise misinterpret a partially
% flushed line (captured right before a timeout kill) as real, if
% incomplete, TISEAN output, producing a silently wrong feature value
% instead of a diagnosable failure. Any other exit code means the TISEAN
% binary itself ran to completion (whether cleanly or with its own
% terse refusal, e.g. false_nearest's exit 54 for "too large" -- see
% NL_FNN.m) and is passed through for the caller to interpret.
%
% The command's output is NOT taken from system()'s own captured output. On
% Linux, under load (e.g. many parallel pool workers on one node), MATLAB's
% system() can return before it has read all of the child's output, and deliver
% the remainder as the start of the NEXT system() call's output in the same
% process. A TISEAN call then returns truncated output (a missing final block,
% which parses without complaint, or an error), and the call after it gets a
% stray numeric tail in front of its own output. The output is therefore
% collected through a file and checked for completeness: see BF_SystemCaptured.m.
%
%---INPUTS:
% tiseanCommand, the TISEAN command-line string to run (as would be
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

persistent timeoutCmd timeoutChecked
if isempty(timeoutChecked)
    % Cache which (if either) timeout utility is on the path -- checked
    % once per MATLAB session, not on every call.
    %
    % The probe runs '<cmd> --version' and inspects ONLY the exit status
    % (0 = present). This deliberately uses no shell redirection or shell
    % builtins so it behaves identically whichever shell MATLAB's system()
    % invokes (it follows $SHELL, which may be csh/tcsh, where sh-style
    % syntax such as 'command -v' or '2>&1' fails).
    [statusTimeout,~]  = system('timeout --version');
    [statusGtimeout,~] = system('gtimeout --version');
    if statusTimeout == 0
        timeoutCmd = 'timeout';
    elseif statusGtimeout == 0
        timeoutCmd = 'gtimeout'; % GNU coreutils' macOS/Homebrew name
    else
        timeoutCmd = ''; % neither available -- degrade to unbounded system()
        warning('BF_TiseanSystem:noTimeoutUtility', ...
            ['No ''timeout'' or ''gtimeout'' utility found on the path: TISEAN ' ...
             'calls will run unbounded, so a hung or runaway TISEAN binary can ' ...
             'stall the computation or exhaust memory. Install GNU coreutils to ' ...
             'restore the wall-clock guard. (This message prints once per session.)']);
    end
    timeoutChecked = true;
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

% A TISEAN call can fail to start, or be killed, for reasons unrelated to the
% data: the shell or the dynamic loader failing transiently under load (status
% 126/127 or a negative status for a command that is installed, or the
% executable, its libraries or its input file briefly unreadable on a busy
% network filesystem), the process being killed by a signal (status 129-192,
% e.g. 137 = SIGKILL from the memory cgroup, 139 = SIGSEGV), a non-zero exit
% that printed nothing at all, or the exit marker (see the header) being absent.
% Such a call is repeated a few times, after a short randomized pause, before it
% is reported as a failure: re-running a deterministic TISEAN binary on the same
% input file either gives the same answer or fixes the glitch, so a retry cannot
% change any value. Statuses the binaries choose themselves (e.g. 54, 50: 50-99,
% see tisean_cec.h) and the wrapper's own timeout (124) are never retried.
maxAttempts = 3;
for attempt = 1:maxAttempts
    [status, res, gotMarker] = BF_SystemCaptured(fullCommand);
    isTransient = ~gotMarker || any(status == [126, 127]) || status < 0 || (status > 128 && status < 193) ...
        || (status ~= 0 && isempty(strtrim(res)) && ~(status >= 50 && status <= 99)); % exit with no message at all
    if ~isTransient || attempt == maxAttempts
        break
    end
    pause(0.2*attempt + 0.3*rand);
end

if clearLdLibraryPath
    setenv('LD_LIBRARY_PATH', oldLdLibraryPath);
end

if ~gotMarker
    error('BF_TiseanSystem:tiseanFailed', ...
        ['TISEAN command did not report an exit status (the shell running it died or could not write its output file; ' ...
         '%d attempt(s)).\nCommand: %s\nOutput: %s'], attempt, tiseanCommand, res);
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
        'TISEAN command could not be completed (%s; exit status %d after %d attempt(s)).\nCommand: %s\nOutput: %s', ...
        reason, status, attempt, tiseanCommand, res);
end

% A call that was still killed by a signal after the retries is passed on to the
% caller exactly as before (same status and output, so no feature value changes),
% but with a warning rather than silently:
if status > 128 && status < 193
    warning('BF_TiseanSystem:tiseanKilled', ...
        'TISEAN command was terminated by a signal (exit status %d, %d attempt(s)): %s', ...
        status, attempt, tiseanCommand);
end

% Belt-and-braces guard against a runaway binary -- one that streams console
% output without terminating (a mode 'lyap_r' and 'false_nearest' can enter
% on pathological input). A working 'timeout' bounds this by wall-clock, and
% every legitimate TISEAN capture here is a few KB at most (numeric results
% mostly go to a '-o' file, not stdout), so an implausibly large 'res' means
% the call misbehaved: treat it as a failure of the same class as a timeout
% rather than handing megabytes of junk to the caller's numeric parser. Also
% covers the degraded no-timeout-utility path above, where the only other
% bound on output volume is available memory.
maxResChars = 50e6; % ~50 MB of text
if numel(res) > maxResChars
    error('BF_TiseanSystem:runawayOutput', ...
        ['TISEAN command produced %d characters of output (limit %d) without ' ...
         'failing cleanly -- treating as a runaway/hung call.\nCommand: %s'], ...
        numel(res), maxResChars, tiseanCommand);
end

end
