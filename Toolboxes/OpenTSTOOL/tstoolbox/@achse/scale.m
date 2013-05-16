function r = scale(a,f)

%tstoolbox/@achse/scale
%   Syntax:
%     * r = scale(a,f)
%
%   Scale achses delta by factor f.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if isa(f, 'double')
	r = a;
	r.delta = f * a.delta;
	r.values = f * a.values;
else
	error('need scalar factor to scale achse')
end
