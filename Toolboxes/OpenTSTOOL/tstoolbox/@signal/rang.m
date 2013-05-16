function rs = rang(s)

%tstoolbox/@signal/rang
%   Syntax:
%     * rs = rang(s)
%
%   Transform scalar time series to rang values.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

error(nargchk(1,1,nargin));

c = rang(s.core); 
rs = signal(c, s);				% special constructor calling syntax for working routines


rs = setyunit(rs, unit);
rs = setyname(rs, 'Rang');
rs = addhistory(rs,  ['Transformed to rang values'] );
rs = addcommandlines(rs, 's = rang(s');
