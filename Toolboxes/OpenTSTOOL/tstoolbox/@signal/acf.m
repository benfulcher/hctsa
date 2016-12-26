function rs = acf(s, len)

%tstoolbox/@signal/acf
%   Syntax:
%     * acf(s, len)
%
%   Input arguments:
%     * len -length of the fft (optional)
%
%   Autocorrelation function for real scalar signals, using fft (of length
%   len). If len is ommited a default value is calculated. The maximum of
%   the calculated length is 128.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);
     
if (ndim(s) > 1) | (~isreal(data(s)))
	help(mfilename)
	return
end

if (nargin < 2) 
	if dlens(s,1) > 256
		len = 128;
	else
		len = nextpow2(dlens(s,1)/4);
	end
end

c = acf(s.core, len);
rs = signal(c, s);	% special constructor calling syntax for working routines
a = getaxis(s, 1); 
dl = delta(a);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit);		% acf values are scalars without unit
rs = addhistory(rs, 'Autocorrelation function');
rs = setlabel(rs, 'Autocorrelation function');
rs = addcommandlines(rs, 's = acf(s', len);
