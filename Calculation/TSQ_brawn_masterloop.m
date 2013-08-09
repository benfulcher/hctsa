function [masteroutput, mastertime] = TSQ_brawn_masterloop(x, y, Mcode, fid, bevocal)

if bevocal
    % Display code name for error checking
    fprintf(fid,'%s...',Mcode);
end

try
	mastertimer = tic;
    if bevocal
        % Any output text is printed to screen
    	masteroutput = BF_pareval(x,y,Mcode,1);
    else
        % Output text stored in T (could log this if you really want to)
        [masteroutput, T] = BF_pareval(x,y,Mcode,0);
    end
	mastertime = toc(mastertimer);
    if bevocal
        fprintf(1,' evaluated (%s).\n',BF_thetime(mastertime))
    end
	% For not-applicable/'real NaN', masteroutput is a NaN, otherwise a
	% structure with components to be called below by pointer operations.
catch emsg
    if bevocal
        fprintf(1,' error (%s).\n',BF_thetime(mastertime))
    end
	fprintf(fid,'---Error evaluating %s\n',Mcode);
    fprintf(fid,'%s\n',emsg.message)
    masteroutput = {}; % Keep empty output
    mastertime = 0; % Set zero calculation time
	% Remains an empty cell entry.
end

end