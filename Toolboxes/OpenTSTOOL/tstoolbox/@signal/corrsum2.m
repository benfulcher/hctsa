function [rs, s] = corrsum(s, npairs, range, past, bins, opt_flag)

%tstoolbox/@signal/corrsum2
%   Syntax:
%     * rs = corrsum2(s, npairs, range, past, bins)
%
%   Input arguments:
%     * npairs - number of pairs per bins
%     * range - maximal relative search radius (relative to attractor
%       size) 0..1
%     * past - number of samples to exclude before and after each
%       reference index
%     * bins - number of bins (optional), defaults to 32
%
%   Compute scaling of correlation sum for time-delay reconstructed
%   timeseries s (Grassberger-Proccacia Algorithm), using fast nearest
%   neighbor search.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,6)

if nargin < 5
	bins  = 32;
end
if nargin < 6
	opt_flag = 1; % to suppress output?
end

[N,dim] = size(data(s));

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

[c,d] = corrsum2(atria, data(s), npairs, range, past, bins, opt_flag);		% call mex-file 'corrsum'

d = log2(d);

rs = signal(core(log2(c(:))), s);	
a = achse(unit, d(1), mean(diff(d)));
a = setname(a, 'ln r');
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed correlation dimension (GPA)']);
rs = addcommandlines(rs, 's = corrsum2(s', npairs, range, past, opt_flag);
rs = setyname(rs, 'ln C(r)');
rs = setlabel(rs, 'Correlation sum');

end
