function [E1, E2] = cao(s, maxdim, tau, NNR, Nref) 

%tstoolbox/@signal/cao
%   Syntax:
%     * [E1, E2] = cao(s, maxdim, tau, NNR, Nref)
%
%   Input arguments:
%     * s - scalar input signal
%     * maxdim - maximal dimension
%     * tau - delay time
%     * NNR - number of nearest neighbor to use
%     * Nref - number of reference points (-1 means: use all points)
%
%   Estimate minimum embedding dimension using Cao's method.
%
%   The second output argument, E2, can be used to distinguish between
%   deterministic and random data.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(5,5);

if (ndim(s) > 1) | (~isreal(data(s)))
	help(mfilename)
	return
end

e = embed(s, maxdim+2, tau);
L = dlens(e, 1);

ref = randref(2, L-2, Nref);
[E, Estar] = cao(data(e), ref, NNR); % calls the mex vesion of cao

c = core(E(2:end)./E(1:end-1));
E1 = signal(c, s);      % special constructor calling syntax for working routines
a = achse(unit, 1, 1);
a = setname(a, 'Dimension (d)');
E1 = setaxis(E1, 1, a);
E1 = setyunit(E1, unit);                % cao values are scalars without unit
E1 = setlabel(E1, 'Minimum embedding dimension using Cao''s method');
E1 = setyname(E1, 'E1(d)');
E1 = addhistory(E1, 'E1 : Minimum embedding dimension estimation using Cao''s method');
E1 = addcommandlines(E1, 's = cao(s', maxdim, tau, NNR, Nref); 

c = core(Estar(2:end)./Estar(1:end-1));
E2 = signal(c, s);      % special constructor calling syntax for working routines
a = achse(unit, 1, 1);
a = setname(a, 'Dimension (d)');
E2 = setaxis(E2, 1, a);
E2 = setyunit(E2, unit);                % cao values are scalars without unit
E2 = addhistory(E2, 'E2 : Minimum embedding dimension estimation using Cao''s method');
E2 = addcommandlines(E2, 's = cao(s', maxdim, tau, NNR, Nref); 
E2 = setlabel(E2, 'Minimum embedding dimension using Cao''s method');
E2 = setyname(E2, 'E2(d)');

end