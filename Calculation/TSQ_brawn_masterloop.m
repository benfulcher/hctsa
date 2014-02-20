% ------------------------------------------------------------------------------
% TSQ_brawn_masterloop
% ------------------------------------------------------------------------------
% 
% Function used in a loop by TSQ_brawn to evaluate a given master function.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

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
        fprintf(1,' error.\n') % ,BF_thetime(mastertime)
    end
	fprintf(fid,'---Error evaluating %s\n',Mcode);
    fprintf(fid,'%s\n',emsg.message)
    masteroutput = {}; % Keep empty output
    mastertime = 0; % Set zero calculation time
	% Remains an empty cell entry.
end

end