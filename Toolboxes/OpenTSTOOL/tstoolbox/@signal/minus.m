function rs = minus(s1,s2)

%tstoolbox/@signal/minus
%   Syntax:
%     * rs=minus(s, offset)
%     * rs=minus(s1,s2)
%
%   Input arguments:
%     * s, s1, s2 - signal object
%     * offset - scalar value
%
%   Calculate difference of signals s1 and s2 or substract a scalar value
%   from s.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(s2, 'signal')
	c = s1.core - s2.core;
%	d = merge(s1.description, s2.description);
	rs = signal(c, s1);
%	rs.description = d;
	rs = addhistory(rs,  ['Subtracted second signal from first signal']);
	rs = addcommandlines(rs, 's = minus(s, s2');
elseif isa(s2, 'double') 
	if size(s2)==[1 1]
		c = core(data(s1)-s2);
		rs = signal(c, s1);
		rs = addhistory(rs,  ['Subtracted ' num2str(s2) ' from signal']);
		rs = addcommandlines(rs, 's = minus(s', s2);
	else
	 	error('second argument must be a signal or a scalar');
	end
else
	error('second argument must be a signal or a scalar');
end

