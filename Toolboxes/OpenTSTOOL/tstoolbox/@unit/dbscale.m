function sc =  dbscale(q)

%tstoolbox/@unit/dbscale
%   returns scaling value when calculating decibel values from data of
%   this unit. dpscale returns either 10 (for power or energy units (e.g.
%   Watt)) or 20 (for all other units (e.g. Volt).
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

sc = q.dBScale;
