function rs = surrogate3(s)

%tstoolbox/@signal/surrogate3
%   Syntax:
%     * rs = surrogate3(s)
%
%   create surrogate data for a scalar time series by permuting samples
%   randomly
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = surrogate3(s.core); 		% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines

rs = addhistory(rs,  ['Surrogated by random sample permutation'] );
rs = addcommandlines(rs, 's = surrogate3(s');

