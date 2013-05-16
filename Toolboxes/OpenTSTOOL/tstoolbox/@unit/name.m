function f =  name(q)

%tstoolbox/@unit/name
%   returns name of unit q.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if isa(q, 'unit')
	f = q.name;
else
	f = '';
end
