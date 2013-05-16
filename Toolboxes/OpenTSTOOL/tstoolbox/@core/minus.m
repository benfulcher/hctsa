function r = minus(c1,c2)

%tstoolbox/@core/minus
%   Syntax:
%     * minus(c1,c2)
%
%   Input Arguments:
%     * c1,c2 - core objects
%
%   subtract c2 from each columns of c1
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if (ndim(c2) == 1) & (ndim(c1)>1)	% subtract c2 from each columns of c1
	dl = dlens(c1);
	r = core(data(c1) - repmat(data(c2), [1 dl(2:end)]));	
else
	r = core(data(c1) - data(c2));
end
