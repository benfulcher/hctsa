function v = spacing(s, dim)

%tstoolbox/@signal/spacing
%     * v = spacing(s) (dim=1)
%     * v = spacing(s, dim)
%
%   return spacing values for xaxis nr. dim
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,2);

dlen = dlens(s);

if nargin == 1
        v = spacing(getaxis(s, 1), dlen(1));
else
        if (dim < 1) | (dim > ndim(s))
	      error(['dim must be between 1 and' num2str(ndim(s))]);
        end
        v = spacing(getaxis(s, dim), dlen(dim));	
end

