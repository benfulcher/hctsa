function rs = ccf(s1, s2, len)

%tstoolbox/@signal/ccf
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);
if (ndim(s1) > 2) | (~isreal(data(s1)))
	help(mfilename)
	return
end

if (nargin < 3) 
	if dlens(s1,1) > 256
		len = 128;
	else
		len = nextpow2(dlens(s1,1)/4);
	end
end

c = xcorr(data(s1), data(s2),'coeff');
rs = signal(c, s1);	% special constructor calling syntax for working routines
a = getaxis(s1, 1); 
dl = delta(a);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = addhistory(rs, 'Cross correlation function');
rs = addcommandlines(rs, 's = ccf(s, s2', len);
