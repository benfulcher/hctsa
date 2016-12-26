function [rs, s] = infodim2(s, n, kmax, past)

%tstoolbox/@signal/infodim2
%   Syntax:
%     * rs = infodim2(s, n, kmax, past)
%
%   Input arguments:
%     * n - number of randomly chosen reference points (n == -1 means :
%       use all points)
%     * kmax - maximal number of neighbors for each reference point
%     * past - number of samples to exclude before and after each
%       reference index
%
%   Compute scaling of moments of the nearest neighbor distances for
%   time-delay reconstructed timeseries s. This can be used to calculate
%   information dimension D1.
%
%   Numerically compute first derivative of logg(k) after k.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,4)


[N,dim] = size(data(s));

ref = randref(1,N, n);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

[nn, dist] = nn_search(data(s), atria, ref, kmax, past);

d = mean(log(dist))';	%'
info = dloggamma(1:kmax); 

rs = signal(core(info-log(N)), s);	
a = achse(unit, d);
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed information dimension']);
rs = addcommandlines(rs, 's = infodim2(s',  n, kmax, past);


function v = dloggamma(k)

% numerically compute first derivative of log(gamma(k)) after k

h = 0.001;
k = [k(:)-h k(:)+h];
v = (gammaln(k(:,2)) - gammaln(k(:,1))) / (2*h);
