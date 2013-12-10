function [opoutput, opquality, optime] = TSQ_brawn_oploop(MasterOutput, MasterCalcTime, MasterLabel, OperationCode, fid)

% Make links to the master
% MasterOutput = MasterOutput{moplink};
% MasterCalcTime = MasterCalcTime(moplink);
% MasterLabel = MasterLabel{moplink};
    
try
    if iscell(MasterOutput) % there was an error evaluating this master operation:
        % fprintf(1,'***Error evaluating master operation %s\n',MasterLabel)
        opoutput = 0; % Output = 0
		opquality = 1; % fatal error QualityCode -- this should not have happened
        optime = NaN; % don't worry about calculation time for errors
        
    elseif ~isstruct(MasterOutput)
        if isnan(MasterOutput); % output from Master was a NaN
    		opoutput = NaN; % all structure elements set to NaN
            optime = NaN; % calculation times also set to NaN
            opquality = NaN; % quality will be set later
        else % A single output -- retrieve it
            opoutput = MasterOutput;
            opquality = 0; % assume it's good -- later special value codes will be assigned where appropriate
            optime = MasterCalcTime;
        end
        
	else % Retrieve the required element from master structure
        [~, thest] = strtok(OperationCode,'.'); thest = thest(2:end); % the field after the '.'
        opoutput = BF_parevalM(MasterOutput,['themasterdat.', thest]);
		opquality = 0; % no evaluation error, quality = 0
        optime = MasterCalcTime;
	end
    
catch emsg
	fprintf(fid,'-----Error linking to master operation %s by %s\n',MasterLabel,OperationCode);
    fprintf(fid,'%s\n',emsg.message)
    keyboard
    opoutput = 0; % Output = 0
	opquality = 1; % fatal error QualityCode -- something of a different error, though...
    optime = NaN; % don't worry about calculation time for errors
end

end