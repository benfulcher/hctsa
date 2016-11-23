function rs = medianfilt(s, len)

%tstoolbox/@signal/medianfilt
%   Syntax:
%     * rs = medianfilt(s, len)
%
%   Moving median filter of width len samples for a scalar time series
%   (len should be odd).
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


narginchk(2,2);

c = medianfilt(s.core, len);
rs = signal(c, s);				% special constructor calling syntax for working routines


a = getaxis(s,1);
a = setfirst(a, ((len+1)/2) * delta(a));

rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Moving average with window of ' num2str(len) ' samples'] );
rs = addcommandlines(rs, 's = medianfilt(s', len);
