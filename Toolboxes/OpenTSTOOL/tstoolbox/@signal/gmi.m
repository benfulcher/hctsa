function rs = gmi(s, D, eps, NNR, len, Nref)

%tstoolbox/@signal/gmi
%   Syntax:
%     * gmi(s, D, eps, NNR, len, Nref)
%
%   Input arguments:
%     * D -
%     * eps -
%     * NNR -
%     * len -
%     * Nref -
%
%   Generalized mutual information function for a scalar time series
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt
narginchk(6,6);


if (ndim(s) > 1) | (~isreal(data(s)))
	help(mfilename)
	return
end


L = dlens(s, 1);
[y,i] = sort(data(s));
ra(i) =  0:(1/L):(1-1/L);		% create rank values

e = embed(signal(ra'), D, 1);		%'

c = core(gmi(data(e), ceil(rand(Nref,1)*(L-len-D-3))+2, eps, NNR, len));
 
rs = signal(c, s);      % special constructor calling syntax for working routines
a = getaxis(s, 1); 
dl = delta(a);
a = setfirst(a, 0);
rs = setaxis(rs, 1, a);
rs = setyunit(rs, unit);                % gmi values are scalars without unit
rs = addhistory(rs, 'Generalized mutual information function');
rs = addcommandlines(rs, 's = gmi(s', D, eps, NNR, len, Nref);
