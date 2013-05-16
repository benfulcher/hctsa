function sspe = mixed (s, dims, future, lags, normhandle, shifts, ...
		       exclude,mode, nn_number,exp)
%
% MIXED mixed states correlation method
%
% 17.6.2002 - Kevin Bube  Drittes Physikal. Institut Uni Goettingen
%
% Copyright 1997-2002 DPI Goettingen
% License http://www.physik3.gwdg.de/tstool/gpl.txt


c = mixed (s.core);
sspe = signal (c, s);

rs = addhistory (rs, ['Single step prediction errors computed by' ...
		    ' mixed state method']);
rs = addcommandlines (rs, ['sspe = mixed (...)']);
