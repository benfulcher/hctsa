function f =  exponents(q)

%tstoolbox/@unit/exponents
%   returns dimension exponents of unit q.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if isa(q, 'unit')
	f = q.exponents;
else
	f = [0 0 0 0 0 0 0 0];
end
