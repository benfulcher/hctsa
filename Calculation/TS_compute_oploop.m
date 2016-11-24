function [opOutput, opQuality, opTime] = TS_compute_oploop(masterOutput, masterCalcTime, masterLabel, operationCode)
% TS_compute_oploop     Links operations to outputs of corresponding master operations,
% including the assignment of error codes.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

try
    if iscell(masterOutput) % there was an error evaluating this master operation:
        % fprintf(1,'***Error evaluating master operation %s\n',masterLabel)
        opOutput = 0; % Output = 0
		opQuality = 1; % fatal error QualityCode -- this should not have happened
        opTime = NaN; % don't worry about calculation time for errors

    elseif ~isstruct(masterOutput)
        if isnan(masterOutput); % output from Master was a NaN
    		opOutput = NaN; % all structure elements set to NaN
            opQuality = NaN; % quality will be set later
            opTime = NaN; % calculation times also set to NaN
        else % A single output -- retrieve it
            if isempty(masterOutput)
                opOutput = 0;
                opQuality = 6;
                opTime = NaN;
            else
                opOutput = masterOutput;
                opQuality = 0; % assume it's good -- later special value codes will be assigned where appropriate
                opTime = masterCalcTime;
            end
        end

	else % Master code returned a structure: retrieve the required element from the master structure:
        % First isolate the field after the '.': theField
        [~, theField] = strtok(operationCode,'.');
        theField = theField(2:end);
        opOutput = masterOutput.(theField);
        if isempty(opOutput) % the field is empty
            opOutput = 0;
            opQuality = 6; % label indicates empty
            opTime = NaN;
        else
    		opQuality = 0; % No evaluation error, quality = 0 means good.
            opTime = masterCalcTime;
        end
	end

catch emsg
    fprintf(1,'-----Error linking to master operation %s by %s\n',masterLabel,operationCode);
    fprintf(1,'%s\n',emsg.message);
    opOutput = 0; % Output = 0
	opQuality = 7; % fatal error QualityCode -- something of a different error, though...
    opTime = NaN; % don't worry about calculation time for errors
end

end
