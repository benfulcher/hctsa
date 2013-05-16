function a = cut(a, start, varargin)

%tstoolbox/@achse/cut
%   Syntax:
%     * a = cut(a, start, stop)
%
%   Cut a part out of achse a, beginning from index start up to index
%   stop. stop is only needed in case of arbitrary spacing. cut ensures
%   the following:
%If values = spacing(achse1, N) and N > n then
%values(n:N) == spacing(cut(achse1, n), N+1-n)
%
%   See also: horzcat
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


switch  resolution(a)
    case 'linear'
		a = setfirst(a, first(a) + (start-1)*delta(a));
    case 'logarithmic'
		a = setfirst(a, first(a) * delta(a) .^ (start-1));	
	case 'arbitrary'
		stop = varargin{1};
		a.values = a.values(start:stop);
end
