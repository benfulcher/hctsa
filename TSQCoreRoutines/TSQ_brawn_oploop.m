function [opoutput, opquality, optime] = TSQ_brawn_oploop(x, y, moplink, Moutput, Mct, Mmlab, parmcodej, fid)

% make links to the master
Moutput = Moutput{moplink};
Mct = Mct(moplink);
Mmlab = Mmlab{moplink};
    
try
    if iscell(Moutput) % there was an error evaluating this master operation:
        % fprintf(1,'***Error evaluating master operation %s\n',Mmlab)
        opoutput = 0; % Output = 0
		opquality = 1; % fatal error QualityCode -- this should not have happened
        optime = NaN; % don't worry about calculation time for errors
    elseif ~isstruct(Moutput)
        if isnan(Moutput); % output from Master was a NaN
    		opoutput = NaN; % all structure elements set to NaN
            optime = NaN; % calculation times also set to NaN
            opquality = NaN; % quality will be set later
        else % A single output -- retrieve it
            opoutput = Moutput;
            opquality = 0; % assume it's good -- later special value codes will be assigned where appropriate
            optime = Mct;
        end
	else % Retrieve required element from master structure
        [~, thest] = strtok(parmcodej,'.'); thest = thest(2:end); % the field after the '.'
        opoutput = BF_parevalM(Moutput,['themasterdat.', thest]);
		opquality = 0; % no evaluation error, quality = 0
        optime = Mct;
	end
catch emsg
	fprintf(fid,'-----Error linking to master operation %s by %s\n',Mmlab,parmcodej);
    fprintf(fid,'%s\n',emsg.message)
    keyboard
    opoutput = 0; % Output = 0
	opquality = 1; % fatal error QualityCode -- something of a different error, though...
    optime = NaN; % don't worry about calculation time for errors
end

end