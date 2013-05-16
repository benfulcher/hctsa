function nd = ndim(c)

%tstoolbox/@core/ndim
%   Syntax:
%     * ndim(c)
%
%   Input Arguments:
%     * c - core object
%
%   return number of dimensions, a scalar value has 0 dimensions
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

nd = length(c.dlens);
if nd == 1
	if c.dlens(1) == 1
		nd = 0;
	end
end
