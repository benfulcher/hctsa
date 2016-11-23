function rs = abs(s)

%tstoolbox/@signal/abs
%   Syntax:
%     * abs(s)
%
%   Take absolut value of all data values of signal s. If sample values
%   are complex, abs(s) returns the complex modulus (magnitude) of each
%   sample.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(1,1);

rs = signal(core(abs(data(s))), s);	
rs = addhistory(rs, ['Absolut values']);
rs = addcommandlines(rs, 's = abs(s');
