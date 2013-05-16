function s = get(l, nr)

%tstoolbox/@list/get
%   Syntax:
%     * s = get(l, nr)
%
%   returns string number nr from list l.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


if (nr >=1) & (nr <= l.len)
	s = l.data{nr};	
else
	s = '';
end
