function d = setyunit(d, u)

%tstoolbox/@description/setyunit
%   Syntax:
%     * d = setyunit(d, unit)
%     * d = setyunit(d, string)
%
%   Set signal's y-unit
%   e.g. d = setyunit(d, 'V')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(u, 'unit')
	d.yunit = u;
elseif isa(u, 'char')
	d.yunit = unit(u);
else
	error('Wrong type of argument(s) given');
end
