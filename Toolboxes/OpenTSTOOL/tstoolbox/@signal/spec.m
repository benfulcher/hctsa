function rs = spec(s)

%tstoolbox/@signal/spec
%   Syntax:
%     * rs = spec(s)
%
%   compute power spectrum for real valued scalar signals. Multivariate
%   signals are accepted but may produce unwanted results as only the
%   spectrum of the first column is returned.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,1);

c = spec(s.core); 	% call real working routine for parent core object
rs = signal(c, s);	% special constructor calling syntax for working routines
   
a = getaxis(s, 1); 
rs = setaxis(rs, 1, achse(unit(a)^(-1),0, samplerate(a)/dlens(s,1)));
rs = addhistory(rs, 'Calculated spectrum (spec)');
rs = setyunit(rs, yunit(s)^2);
rs = addcommandlines(rs, 's = spec(s');
