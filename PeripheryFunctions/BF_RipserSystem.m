function [status, res] = BF_RipserSystem(ripserCommand, timeoutSecs)
% BF_RipserSystem   Run a ripser command-line call with a bounded timeout.
%
% Sibling of BF_TiseanSystem.m for the ripser persistent-homology binary
% (Toolboxes/ripser), with the same timeout-wrapping and exit-code
% classification logic -- kept as a separate function rather than a shared
% generic wrapper so that this addition doesn't touch BF_TiseanSystem.m, an
% already-tested shared dependency of ~16 existing operations.
%
% A call that fails this way -- status 124 (this wrapper's own timeout), or
% 126/127 (the shell's reserved "found but not executable"/"command not
% found" codes, i.e. ripser never ran at all) -- errors immediately rather
% than returning res to the caller, for the same reason as BF_TiseanSystem:
% a caller parsing persistence intervals out of res could otherwise
% misinterpret partially flushed output as a complete (if small) diagram.
% Any other exit code means ripser itself ran to completion and is passed
% through as before for the caller to interpret.
%
%---INPUTS:
% ripserCommand, the ripser command-line string to run (as previously passed
%                 directly to system())
% timeoutSecs [opt], wall-clock time limit in seconds (default: 600)
%
%---OUTPUTS:
% status, the exit status, as returned by system() for the ripser binary's
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
    % Cache which (if either) timeout utility is on the path -- checked once
    % per MATLAB session, not on every call (see BF_TiseanSystem.m for why
    % 'timeout' vs 'gtimeout' matters on macOS):
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
    [status, res] = system(ripserCommand);
else
    [status, res] = system(sprintf('%s %g %s', timeoutCmd, timeoutSecs, ripserCommand));
end

if any(status == [124, 126, 127])
    switch status
        case 124
            reason = sprintf('timed out after %g seconds', timeoutSecs);
        case 126
            reason = 'found but not executable';
        case 127
            reason = 'not found on the system path -- is ripser installed and compiled? (see install.m / compile_ripser.m)';
    end
    error('BF_RipserSystem:ripserFailed', ...
        'ripser command could not be completed (%s).\nCommand: %s\nOutput: %s', ...
        reason, ripserCommand, res);
end

end
