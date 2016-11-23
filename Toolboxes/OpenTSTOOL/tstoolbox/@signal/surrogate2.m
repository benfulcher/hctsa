function rs = surrogate2(s)

%tstoolbox/@signal/surrogate2
%   Syntax:
%     * rs = surrogate2(s)
%
%   create surrogate data for a scalar time series
%   see : James Theiler et al.'Using Surrogate Data to Detect Nonlinearity
%   in Time Series', APPENDIX : ALGORITHM II
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = surrogate2(s.core); 		% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines

rs = addhistory(rs,  ['Surrogated with Theiler Algorithm II'] );
rs = addcommandlines(rs, 's = surrogate2(s');
