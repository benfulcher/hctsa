function f =  quantity(q, which)

%tstoolbox/@unit/quantity
%   returns quantity name of unit q. If argument which is omitted, the
%   english quantity name will be returned.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if isa(q, 'unit')
	if nargin > 1
		if strcmp(which, 'ger')
			f = q.quantity.ger;
		else
			f = q.quantity.eng;
		end
	else
		f = q.quantity.eng;	
	end
else
	f = '';
end
