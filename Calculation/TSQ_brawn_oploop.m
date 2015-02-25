% ------------------------------------------------------------------------------
% TSQ_brawn_oploop
% ------------------------------------------------------------------------------
% 
% Function that links operations to outputs of corresponding master operations,
% including the assignment of error codes.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function [opOutput, opQuality, opTime] = TSQ_brawn_oploop(masterOutput, masterCalcTime, masterLabel, operationCode, fid)
    
try
    if iscell(masterOutput) % there was an error evaluating this master operation:
        % fprintf(1,'***Error evaluating master operation %s\n',masterLabel)
        opOutput = 0; % Output = 0
		opQuality = 1; % fatal error QualityCode -- this should not have happened
        opTime = NaN; % don't worry about calculation time for errors
        
    elseif ~isstruct(masterOutput)
        if isnan(masterOutput); % output from Master was a NaN
    		opOutput = NaN; % all structure elements set to NaN
            opTime = NaN; % calculation times also set to NaN
            opQuality = NaN; % quality will be set later
        else % A single output -- retrieve it
            opOutput = masterOutput;
            opQuality = 0; % assume it's good -- later special value codes will be assigned where appropriate
            opTime = masterCalcTime;
        end
        
	else % Retrieve the required element from master structure
        [~, thest] = strtok(operationCode,'.'); thest = thest(2:end); % the field after the '.'
        opOutput = BF_parevalM(masterOutput,['themasterdat.', thest]);
		opQuality = 0; % No evaluation error, quality = 0 means good.
        opTime = masterCalcTime;
	end
    
catch emsg
	fprintf(fid,'-----Error linking to master operation %s by %s\n',masterLabel,operationCode);
    fprintf(fid,'%s\n',emsg.message)
    keyboard
    opOutput = 0; % Output = 0
	opQuality = 1; % fatal error QualityCode -- something of a different error, though...
    opTime = NaN; % don't worry about calculation time for errors
end

end