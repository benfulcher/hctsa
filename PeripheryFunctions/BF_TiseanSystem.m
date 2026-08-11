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
%---INPUTS:
% tiseanCommand, the TISEAN command-line string to run (as previously
%                 passed directly to system())
% timeoutSecs [opt], wall-clock time limit in seconds (default: 120)
%
%---OUTPUTS:
% status, the exit status, as returned by system(). If the call is actually
%          killed by the timeout (rather than exiting under its own steam),
%          this is 124 -- GNU timeout's own reserved code for that case, so
%          existing/future status checks can test for it directly.
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
    timeoutSecs = 120;
end

if isempty(timeoutCmd)
    [status, res] = system(tiseanCommand);
else
    [status, res] = system(sprintf('%s %g %s', timeoutCmd, timeoutSecs, tiseanCommand));
end

end
