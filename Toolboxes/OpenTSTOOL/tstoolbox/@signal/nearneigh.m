function [rs, s] = nearneigh(s, n, past)

%tstoolbox/@signal/nearneigh
%   Syntax:
%     * rs = nearneigh(s, n) => past=1
%     * rs = nearneigh(s, n, past)
%
%   Input arguments:
%     * n - number of nearest neighbour to compute
%     * past - a nearest neighbour is only valid if it is as least past
%       timesteps away from the reference point. past = 1 means: use all
%       points but ref_point itself
%
%   n nearest neighbour algorithm. Find n nearest neighbours (in order of
%   increasing distances) to each point in signal s uses accelerated
%   searching, distances are calculated with euclidian norm.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,3);     

if nargin < 3
	past = 1;
end

if ndim(s) ~= 2
	error('need a two-dimensional signal to compute nearest neighbour statistic');
end


past = past;	
points = data(s);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(points, 'euclidian');
	s = setoptparams(s, 1, atria);
end

nn = nn_search(data(s), atria, 1:dlens(s, 1), n, past); 

c = core(nn);
rs = signal(c, s);	% special constructor calling syntax for working routines
rs = setaxis(rs, 2, achse);
rs = setplothint(rs, 'multipoints');
rs = addhistory(rs, {['Calculated indices of the ' num2str(n) ' nearest neighbours, excluding ' num2str(past+1) ' samples before and after the actual index']});
rs = addcommandlines(rs, 'nearneigh(s', n, past);

