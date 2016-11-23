function [xpos, u] = firstmax(s)

%tstoolbox/@signal/firstmax
%   Syntax:
%     * [xpos, unit] = firstmax(s)
%
%   Give information about first local maximum of scalar signal s.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(1,1);

if ndim(s) ~= 1
	help(mfilename)
	return
end

d = data(s);
N = dlens(s,1);

c1 = d(1:N-2);
c2 = d(2:N-1);
c3 = d(3:N);

m = find((c1<c2) .* (c3<c2));
if isempty(m)
	m = find((c1<=c2) .* (c3<c2));
end
if isempty(m)
	m = find((c1<c2) .* (c3<=c2));
end

if isempty(m)
	error('No local minimum found')
else
	A = getaxis(s, 1);
	xpos = first(A) + m(1)*delta(A);
	u = unit(A);
end

