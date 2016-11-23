function s = shift(s, distance, dim)

%tstoolbox/@signal/shift
%   Syntax:
%     * s = shift(s, distance) (dim=1)
%     * s = shift(s, distance, dim)
%
%   shift signal on axis No. dim by distance (measured in the unit of the
%   axis) to the right
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	dim = 1;
end

a = getaxis(s, dim);

switch  resolution(a)
	case 'linear'
		a = setfirst(a, first(a) + distance);
		s = setaxis(s, dim, a);
	otherwise
		error('Data values are not sampled equidistantly');
end

s = addcommandlines(s, 's = shift(s', distance, dim);


