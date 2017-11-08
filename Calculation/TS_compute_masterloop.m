function [masterOutput, masterTime] = TS_compute_masterloop(x, y, masterCode, masterID, numMasterOps, beVocal, theTsID, iterNum)
% TS_compute_masterloop     Used in a loop by TS_compute to evaluate a given master function.

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems (2017).
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

if beVocal
    % Display code name for error checking
    fprintf(1,'[ts_id = %u, mop_id = %u (%u/%u)] %s...', theTsID, masterID, iterNum, numMasterOps, masterCode);
end

try
	masterTimer = tic;
    if beVocal
        % Any output text is printed to screen
    	masterOutput = BF_pareval(x,y,masterCode,1);
    else
        % Output text stored in second output (could log this if you really want to)
        masterOutput = BF_pareval(x,y,masterCode,0);
    end
	masterTime = toc(masterTimer);
    if beVocal
        fprintf(1,' evaluated (%s).\n',BF_thetime(masterTime));
    end
	% For not-applicable/'real NaN', masterOutput is a NaN, otherwise a
	% structure with components to be called below by pointer operations.

catch emsg
    if beVocal
        fprintf(1,' error.\n'); % ,BF_thetime(masterTime)
    end
	fprintf(1,'---Error evaluating %s:\n%s\n',masterCode,emsg.message);
    masterOutput = {}; % Keep empty output
    masterTime = 0; % Set zero calculation time
	% Remains an empty cell entry.
end

end
