function newsig = mixemb (sig, dim1, dim2, delay1, delay2, ...
			  shifts)

% MIXEMB  wrapper function for mixembed
%
% 10.6.2002 - Kevin Bube  Drittes Physikal. Institut Uni Goettingen
%
% Copyright 1997-2002 DPI Goettingen
% License http://www.physik3.gwdg.de/tstool/gpl.txt

tmp    = mixembed (data (sig), [dim1 dim2], [delay1 delay2], shifts);
newsig = signal (tmp);
