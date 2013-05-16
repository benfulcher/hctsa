function rs = int(s)

%tstoolbox/@signal/int
%   Syntax:
%     * int(s)
%
%   Numerical integration along dimension 1 signal s has to be sampled
%   equidistantly.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if nargin < 2
	nth = 1;
end

a = getaxis(s,1);
switch  resolution(a)
	case 'linear'
		c = int(s.core, delta(a)); 		% call real working routine for parent core object
		rs = signal(c, s);				% special constructor calling syntax for working routines
		rs = setyunit(rs, yunit(s)*unit(a));
		a = setfirst(a, first(a) - delta(a)/2);
		rs = setaxis(rs,  1, a);
		rs = addhistory(rs,  ['Numerical integration along dimension 1'] );
		rs = addcommandlines(rs, 's = int(s');
	otherwise
		error('data are not sampled equidistantly');
end


