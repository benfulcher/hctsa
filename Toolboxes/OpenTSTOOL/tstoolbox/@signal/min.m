function [ymin, yu, xpos, xu] = min(s)
 
%tstoolbox/@signal/min
%   Syntax:
%     * [minimum, yunit, xpos, xunit] = min(s)
%
%   Give information about minimum of scalar signal s.
%
%   Example:
%disp('minimum of signal : ')
%disp(['y = ' num2str(m) ' ' label(yunit(s))]);
%disp(['x = ' num2str(xpos) ' ' label(a)]);
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

if ndim(s) ~= 1
	help(mfilename)
	return
end

[m, ind] = min(data(s));
a = getaxis(s, 1);
xpos = first(a) + delta(a) * (ind-1);	% determine position of minimum
% disp('minimum of signal : ')
% disp(['y = ' num2str(m) ' ' label(yunit(s))]);
% disp(['x = ' num2str(xpos) ' ' label(a)]);
xu = unit(a);
yu = yunit(s);
ymin = m;
