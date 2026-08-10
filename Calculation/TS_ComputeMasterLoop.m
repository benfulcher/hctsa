function [masterOutput, masterTime] = TS_ComputeMasterLoop(x,x_z,masterFn,masterCode,masterID,numMasterOps,howVocal,theTsID,iterNum,beVerbose)
% TS_ComputeMasterLoop     Used in a loop by TS_Compute to evaluate a given master function.
%
% masterFn is a function handle, precompiled once from masterCode (see
% TS_CalculateFeatureVector/TS_Compute), taking (x,x_z) as inputs. masterCode
% is kept only for display/error messages below.
%
% beVerbose, [default: false] governs any text (fprintf/disp/warning) that
% masterFn itself prints while it runs: shown, attributed to this master
% operation, when true; silently discarded when false. Independent of
% howVocal, which controls only this function's own progress-reporting text.

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

if nargin < 10 || isempty(beVerbose)
    beVerbose = false;
end

beVocal = strcmp(howVocal,'full') || beVerbose;

if beVocal
    % Display code name for error checking
    fprintf(1,'[TimeSeries_ID = %u, MasterOperation_ID = %u (%u/%u)] %s...',...
                theTsID, masterID, iterNum, numMasterOps, masterCode);
end

try
	masterTimer = tic;
    % Call the precompiled function handle via evalc, which lets us capture (and
    % then either show or discard) any fprintf/disp/warning text masterFn prints
    % internally -- a single point of control for output from arbitrary Operations
    % code, rather than relying on every Operation to gate its own messages.
    % Function handles (rather than eval-ing masterCode as a string) keep this
    % parfor-safe; evalc of the call itself is a plain function call from parfor's
    % point of view and carries negligible overhead (~10us, measured) next to
    % typical feature computation times.
    prevWarnState = warning('off','backtrace'); % keep captured warnings free of clickable stack-trace clutter
    capturedText = evalc('masterOutput = masterFn(x,x_z);');
    warning(prevWarnState);
	masterTime = toc(masterTimer);
    if beVocal
        fprintf(1,' evaluated (%s).\n',BF_TheTime(masterTime));
    end
    if beVerbose && ~isempty(capturedText)
        fprintf(1,'%s',capturedText);
    end
	% For not-applicable/'real NaN', masterOutput is a NaN, otherwise a
	% structure with components to be called below by pointer operations.

catch emsg
    warning(prevWarnState); % restore even though evalc threw mid-capture
    if beVocal
        fprintf(1,' error.\n'); % ,BF_TheTime(masterTime)
    end
	fprintf(1,'---Error evaluating %s:\n%s.\n',masterCode,emsg.message);
    masterOutput = {}; % Keep empty output
    masterTime = 0; % Set zero calculation time
	% Remains an empty cell entry.
end

end
