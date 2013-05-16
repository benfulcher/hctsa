function d = double(u)

%tstoolbox/@unit/double
%   gives a row vector which's first element contains the unit's factor
%   and the remaining elements contain the exponents of the SI basic
%   units.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
 
d = [u.factor u.exponents];
