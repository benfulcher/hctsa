function d = setname(d, n)

%tstoolbox/@description/setname
%   Syntax:
%     * d = setname(d, name)
%
%   the name field of a descriptiom is used when the signal is loaded from
%   file, it will not be continued through several processing steps
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

if isa(n, 'char')
	d.name = n;
else
    error('Wrong type of argument(s) given');
end
