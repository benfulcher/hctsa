function r = eq(a,b)

%tstoolbox/@achse/eq
%   Test if achse a and achse b are equal.first is not (!) taken into
%   account for this test.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

r = ((a.unit == b.unit) & strcmp(a.resolution,b.resolution) & (a.delta == b.delta));
