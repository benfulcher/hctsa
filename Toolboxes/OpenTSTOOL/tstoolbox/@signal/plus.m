function rs = plus(s1,s2)

%tstoolbox/@signal/plus
%   Syntax:
%     * rs=plus(s, offset)
%     * rs=plus(s1, s2)
%
%   Add two signals s1 and s2 or add a scalar value offset to s.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(s2, 'signal')
	c = s1.core + s2.core;
%	d = merge(s1.description, s2.description);
	rs = signal(c, s1);
%	rs.description = s1.description;
	rs = addhistory(rs,  ['Added both signals signal']);
	rs = addcommandlines(rs, 's = plus(s, s2');
elseif isa(s2, 'double') 
	if size(s2)==[1 1]
		c = core(data(s1)-s2);
		rs = signal(c, s1);
		rs = addhistory(rs,  ['Added ' num2str(s2) ' to signal']);
		rs = addcommandlines(rs, 's = plus(s', s2);
	else
	 	error('second argument must be a signal or a scalar');
	end
else
	error('second argument must be a signal or a scalar');
end

