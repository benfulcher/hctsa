function [masteroutput, mastertime] = TSQ_brawn_masterloop(x, y, mindex, Mmcode, fid, bevocal)
	if bevocal
        fprintf(fid,'%s...',Mmcode{mindex});
    end % display for error checking
	try
		mastertimer = tic;
		masteroutput = pareval(x,y,Mmcode{mindex});
		mastertime = toc(mastertimer);
        if bevocal
            fprintf(1,' evaluated.\n')
        end
		% for not-applicable/'real NaN', masteroutput is a NaN, otherwise a
		% structure with components to be called below by pointer operations.
    catch emsg
        if bevocal
            fprintf(1,' error.\n')
        end
		fprintf(fid,'---Error evaluating master operation: %s\n',Mmcode{mindex});
        fprintf(fid,'%s\n',emsg.message)
        masteroutput = {}; % keep empty output
        mastertime = 0; % set zero calculation time
		% remains an empty cell entry.
	end

end