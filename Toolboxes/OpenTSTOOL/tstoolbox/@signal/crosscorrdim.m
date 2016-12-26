function rs = crosscorrdim(s, s2, n, range, past, bins)

%tstoolbox/@signal/crosscorrdim
%   Syntax:
%     * rs = crosscorrdim(s, s2, n, range, past, bins)
%
%   Input arguments:
%     * n - number of randomly chosen reference points (n == -1 means :
%       use all points)
%     * range - maximal relative search radius (relative to size of data
%       set s2) 0..1
%     * past - number of samples to exclude before and after each
%       reference index
%     * bins - number of bins (optional)
%
%   Compute scaling of cross-correlation sum for time-delay reconstructed
%   timeseries s against signal s2 (with same dimension as s), using fast
%   nearest neighbor search. Reference points are taken out of signal s,
%   while neigbors are searched in s2. The default number of bins is 32.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(5,6)

if nargin < 6
	bins  = 32;
end

points = data(s);
[N,dim] = size(points);

points2 = data(s2);
[N2,dim2] = size(points2);

ref = randref(1,N, n);

[c,d] = crosscorrsum(points, points2, ref, range, past, bins);	

d = log2(d);

a = achse(unit, d(1), mean(diff(d)));
a = setname(a, 'ln r');
rs = signal(core(log2(c(:))), s);
newdes = merge(s.description, s2.description);
rs.description = newdes;

rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed cross correlation sum']);
rs = addcommandlines(rs, 's = crosscorrdim2(s,s2', n, range, past);
rs = setyname(rs, 'ln C(r)');
rs = setlabel(rs, 'Cross-correlation sum');
