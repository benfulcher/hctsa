function [rs, s] = localdensity(s, n, past)

%tstoolbox/@signal/localdensity
%   Syntax:
%     * rs = localdensity(s, n, past)
%
%   Input arguments:
%     * n - number of nearest neighbour to compute
%     * past - a nearest neighbour is only valid if it is as least past
%       timesteps away from the reference point past = 1 means: use all
%       points but ref_point itself
%
%   Uses accelerated searching, distances are calculated with euclidian
%   norm.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);

if nargin < 3
	past = 1;
end

if ndim(s) ~= 2
	error('need a two-dimensional signal to compute nearest neighbour statistic');
end

N = dlens(s, 1);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

[nn, dists] = nn_search(data(s), atria, 1:N, n, past); 

c = core(sum(dists,2));
rs = signal(c, s);	% special constructor calling syntax for working routines
rs = setaxis(rs, 2, achse);
rs = addhistory(rs, {['Computed reciprocal local density along trajectory']});
rs = addcommandlines(rs, 'localdensity(s', n, past);

