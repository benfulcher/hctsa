function [rs, s] = return_time(s, nnr, maxT, past, N)

%tstoolbox/@signal/return_time
%   Syntax:
%     * rs = return_time(s, nnr, maxT) => past=1
%     * rs = return_time(s, nnr, maxT, past)
%     * rs = return_time(s, nnr, maxT, past, N)
%
%   Input arguments:
%     * nnr - number of nearest neighbors
%     * maxT - maximal return time to consider
%     * past - a nearest neighbor is only valid if it is as least past
%       timesteps away from the reference point past = 1 means: use all
%       points but tt ref_point itself
%     * N - number of reference indices
%
%   Compute histogram of return times.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(3,5);

if nargin < 4
	past = 1;
end
if nargin < 5
	N = ceil(dlens(s,1)/2);
end

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

hist = return_time(atria, data(s), randref(1,dlens(s,1), N), nnr, maxT, past);

c = core(hist);
rs = signal(c, s);	% special constructor calling syntax for working routines
a = getaxis(rs, 1);
a = setfirst(a, delta(a));
rs = setaxis(rs, 1, a);
rs = addhistory(rs, {['Computed return times']});
rs = addcommandlines(rs, 's = return_time(s', nnr, maxT, past, N);

