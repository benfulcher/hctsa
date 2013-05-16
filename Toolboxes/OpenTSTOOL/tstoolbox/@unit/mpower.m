function r = mpower(u,p)

%tstoolbox/@unit/mpower
%   Syntax:
%     * mpower(u,p)
%
%   take unit u to power p, p must be a scalar.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

if ~isa(p, 'double')
	help(mfilename)
	return
end

if (u.exponents == [0 0 0 0 0 0 0 0]) & (~isempty(u.label))
	r = unit; 
	r.factor = (u.factor)^(p);
	r.label = [u.label '^(' num2str(p) ')']; 
else
	r.exponents = u.exponents * p;
	r.factor = (u.factor)^(p);
	r = unit([(u.factor)^(p) (u.exponents * p)]);
end
