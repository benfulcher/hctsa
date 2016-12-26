function rs = trend(s, len)

%tstoolbox/@signal/trend
%   Syntax:
%     * rs = trend(s, len)
%
%   trend correction
%   calculate moving average of width len (samples) for a scalar time
%   series (len should be odd) and remove the result from the input signal
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,2);

c1 = s.core;
c2 = movav(s.core, len);
N1 = dlens(c1, 1);
N2 = dlens(c2, 1);
shft = floor((N1-N2)/2);
dat = data(c1);
c3 = core(dat(shft:shft+N2-1)); 
c1 = c3 - c2;
rs = signal(c1, s);				% special constructor calling syntax for working routines

a = getaxis(s,1);
a = setfirst(a, (shft-1) * delta(a));
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Removed trend'] );
rs = addcommandlines(rs, 's = trend(s', len);
