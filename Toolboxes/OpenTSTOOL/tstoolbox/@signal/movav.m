function rs = movav(s, len, windowtype)

%tstoolbox/@signal/movav
%   Syntax:
%     * rs = movav(s, len, windowtype)
%     * rs = movav(s, len)
%
%   Moving average of width len (samples) along first dimension.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	c = core(movav(data(s), len));
	rs = signal(c, s);				% special constructor calling syntax for working routines
	rs = addhistory(rs,  ['Moving average with window of ' num2str(len) ' samples'] );
	rs = addcommandlines(rs, 's = movav(s', len);
else
	c = core(movav(data(s), window(len, windowtype)));
	rs = signal(c, s);				% special constructor calling syntax for working routines
	rs = addhistory(rs,  ['Moving average with window of ' num2str(len) ' samples'] );
	rs = addcommandlines(rs, 's = movav(s', len, windowtype);
end

