function [masteroutput, mastertime] = TSQ_brawn_masterloop(x, y, Mcode, fid, bevocal)

if bevocal
    % display code name for error checking
    fprintf(fid,'%s...',Mcode);
end

try
	mastertimer = tic;
	masteroutput = pareval(x,y,Mcode);
	mastertime = toc(mastertimer);
    if bevocal
        fprintf(1,' evaluated (%s).\n',benrighttime(mastertime))
    end
	% for not-applicable/'real NaN', masteroutput is a NaN, otherwise a
	% structure with components to be called below by pointer operations.
catch emsg
    if bevocal
        fprintf(1,' error.\n')
    end
	fprintf(fid,'---Error evaluating master operation: %s\n',Mcode);
    fprintf(fid,'%s\n',emsg.message)
    masteroutput = {}; % keep empty output
    mastertime = 0; % set zero calculation time
	% remains an empty cell entry.
end

end