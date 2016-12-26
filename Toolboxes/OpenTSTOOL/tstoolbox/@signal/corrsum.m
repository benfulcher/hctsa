function [rs, s] = corrsum(s, n, range, past, bins, opt_flag)

%tstoolbox/@signal/corrsum
%   Syntax:
%     * rs = corrsum(s, n, range, past, bins)
%
%   Input arguments:
%     * n - number of randomly chosen reference points (n == -1 means: use
%       all points)
%     * range - maximal relative search radius (relative to attractor
%       size) 0..1
%     * past - number of samples to exclude before and after each
%       reference index
%     * bins - number of bins (optional)
%
%   Compute scaling of correlation sum for time-delay reconstructed
%   timeseries s (Grassberger-Proccacia Algorithm), using fast nearest
%   neighbor search. Default number of bins is 20.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,6)

if nargin < 5
	bins  = 32;
end
if nargin < 6
	opt_flag = 0;
end

[N,dim] = size(data(s));

ref = randref(1,N, n);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

[c,d] = corrsum(atria, data(s), ref, range, past, bins, opt_flag);		% call mex-file 'corrsum'

d = log2(d);

rs = signal(core(log2(c(:))), s);	
a = achse(unit, d(1), mean(diff(d)));
a = setname(a, 'ln r');
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed correlation dimension (GPA)']);
rs = addcommandlines(rs, 's = corrsum(s', n, range, past, bins, opt_flag);
rs = setyname(rs, 'ln C(r)');
rs = setlabel(rs, 'Correlation sum');
