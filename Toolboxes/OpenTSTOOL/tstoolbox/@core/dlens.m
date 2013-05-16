function d = dlens(c, nr)

%tstoolbox/@core/dlens
%   Syntax:
%     * d=dlens(c, nr)
%
%   Input Arguments:
%     * c - core object
%
%   returns sizes of dimensions (same as function 'size' under matlab)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if nargin == 1
	d = c.dlens;
else
	if nr > ndim(c)
		d = 1;
	else
		d = c.dlens(nr);
	end
end
