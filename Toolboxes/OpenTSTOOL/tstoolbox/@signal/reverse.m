function rs = reverse(s)

%tstoolbox/@signal/reverse
%   Syntax:
%     * rs=reverse(s)
%
%   Reverse signal along dimension 1.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

dlen = dlens(s);
points = data(s);

if ndim(s) > 2
	points = reshape(points, [dlen(1) prod(dlen(2:end))]); 
	points = flipud(points);
	points = reshape(points, [dlen(1) dlen(2:end)]); 
else
	points = flipud(points);
end

c = core(points);
rs = signal(c,s);

rs = addhistory(rs, ['Reversed along dimension 1']);
rs = addcommandlines(rs, 's = reverse(s');
