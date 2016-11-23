function rs = surrogate1(s)

%tstoolbox/@signal/surrogate1
%   Syntax:
%     * rs = surrogate1(s)
%
%   create surrogate data for a scalar time series by randomizing phases
%   of fourier spectrum
%   see : James Theiler et al.'Using Surrogate Data to Detect Nonlinearity
%   in Time Series', APPENDIX : ALGORITHM I
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = surrogate1(s.core); 		% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines

rs = addhistory(rs,  ['Surrogated by phase randomization (Theiler Algorithm I)'] );
rs = addcommandlines(rs, 's = surrogate1(s');

