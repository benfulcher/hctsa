function [status, res, gotMarker] = BF_SystemCaptured(command)
% BF_SystemCaptured   Run a shell command, capturing its output through a file
%
% A replacement for [status,res] = system(command) for commands whose output is
% parsed (TISEAN, ripser). On Linux, under load (e.g. many parallel pool workers
% on one node), MATLAB's system() can return before it has read all of the
% child's output, and hand the remainder to the NEXT system() call in the same
% process, as the start of that call's output. A parser then gets a
% truncated result (a missing final block can parse without complaint) or one
% with a stray numeric tail in front.
%
% Here the command runs as
%     sh -c '( <command> 2>&1; printf "\nHCTSA_EXIT=%d\n" $? ) | cat > file'
% so the output goes to a private temporary file, and the shell appends the
% command's own exit status as a last marker line. The command still writes into
% a pipe, exactly as it would under system() (this matters: TISEAN's Fortran
% programs order their stderr banner and their stdout data differently for a pipe
% and for a file, and callers parse the pipe ordering), but `cat` drains that pipe
% to EOF, so the file is complete once the command has exited. system()'s own
% captured stdout is then empty, so nothing can be left behind for a later call.
% Using 'sh -c' also makes the redirection syntax independent of the user's login
% shell (csh/tcsh do not understand 2>&1).
%
%---INPUTS:
% command, the command line to run; it must itself be valid sh syntax
%
%---OUTPUTS:
% status, the command's exit status (as the shell saw it)
% res, the command's stdout and stderr, exactly as it wrote them
% gotMarker, false if the exit-status marker was not found in the output file,
%            i.e. the shell running the command died or could not write the file
%            (status is then system()'s own, and res is whatever was captured)

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

marker = 'HCTSA_EXIT=';
if ispc
    % no sh: use system()'s own capture
    [status, res] = system(command);
    gotMarker = true;
    return
end
outFile = [tempname('/tmp') '.out'];
innerCommand = strrep(command, '''', '''\'''''); % make single quotes safe inside sh -c '...'
shellCommand = sprintf('sh -c ''( %s 2>&1; printf "\\n%s%%d\\n" $? ) | cat > %s''', innerCommand, marker, outFile);
[sysStatus, sysRes] = system(shellCommand);
txt = '';
if isfile(outFile)
    try
        txt = fileread(outFile);
    catch
    end
    try
        delete(outFile);
    catch
    end
end
% The marker (written last, after the command exited) is the final line:
pos = strfind(txt, [newline marker]);
gotMarker = false;
status = sysStatus;
res = txt;
if ~isempty(pos)
    pos = pos(end);
    val = sscanf(txt(pos + 1 + numel(marker):end), '%d', 1);
    if ~isempty(val)
        gotMarker = true;
        status = val;
        res = txt(1:pos - 1);
    end
end
if ~gotMarker
    status = sysStatus;
    res = [txt sysRes];
end

end
