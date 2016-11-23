function d = setyname(d, n)

%tstoolbox/@description/setyname
%   Syntax:
%     * d = setyunit(d, string)
%
%   Set signal's y-name
%   e.g. d = setyunit(d, 'V')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(n, 'char')
	d.yname = n;
else
	error('Wrong type of argument(s) given');
end
