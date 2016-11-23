function rs = power(s)

%tstoolbox/@signal/power
%   Syntax:
%     * power(s)
%
%   Calculate squared magnitude of each sample.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

x = data(s);

rs = signal(core(x .* conj(x)), s);
rs = setyunit(rs, yunit(s)^2);	
rs = addhistory(rs, ['Squared magnitude']);
rs = addcommandlines(rs, 's = power(s');

