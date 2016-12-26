function [ymax, yu, xpos, xu] = max(s)

%tstoolbox/@signal/max
%   Syntax:
%     * [maximum, yunit, xpos, xunit] = max(s)
%
%   Give information about maximum of scalar signal s.
%
%   Example:
%disp('maximum of signal : ')
%disp(['y = ' num2str(m) ' ' label(yunit(s))]);
%disp(['x = ' num2str(xpos) ' ' label(a)]);
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


narginchk(1,1);

if ndim(s) ~= 1
	help(mfilename)
	return
end

[m, ind] = max(data(s));
a = getaxis(s, 1);
xpos = first(a) + delta(a) * (ind-1);	% determine position of maximum
% disp('maximum of signal : ')
% disp(['y = ' num2str(m) ' ' label(yunit(s))]);
% disp(['x = ' num2str(xpos) ' ' label(a)]);
xu = unit(a);
yu = yunit(s);
ymax = m;
