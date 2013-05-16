function s = setunit(s, dim, u)

%tstoolbox/@signal/setunit
%   Syntax:
%     * s = setunir(s, dim, u)
%
%   Change unit of one of the current xaxes.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(3,3,nargin));

a = getaxis(s,dim);
a = setunit(a, u);
s = setaxis(s, dim, a);

