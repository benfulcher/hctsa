function rs = db(s, dbmin)

%tstoolbox/@signal/db
%   Syntax:
%     * db(s, dbmin)
%
%   Compute decibel values of signal relative to a reference value that is
%   determined by the signal's yunit values below dbmin are set to dbmin.
%   If dbmin is ommited it is set to -120.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(1,2);
if nargin < 2, dbmin = -120; end

yu = yunit(s);

ref = dbref(yu);
scf = dbscale(yu);

c = db(s.core, ref, scf, dbmin); 	% call real working routine for parent core object
rs = signal(c, s);				% special constructor calling syntax for working routines

rs = setyunit(rs, unit(['dB' label(yu)]));
rs = addhistory(rs,  ['Calculated decibel values (stretch=' num2str(scf) ',ref=' num2str(ref) ')'] );
rs = addcommandlines(rs, 's = db(s', dbmin);
