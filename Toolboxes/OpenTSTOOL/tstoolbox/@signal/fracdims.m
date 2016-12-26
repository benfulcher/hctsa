function [rs, s] = fracdims(s, kmin, kmax, Nref, gstart, gend, past, steps)

%tstoolbox/@signal/fracdims
%   Syntax:
%     * rs = fracdims(s, kmin, kmax, Nref, gstart, gend, past, steps)
%     * rs = fracdims(s, kmin, kmax, Nref, gstart, gend, past)
%     * rs = fracdims(s, kmin, kmax, Nref, gstart, gend)
%
%   Input arguments:
%     * kmin - minimal number of neighbors for each reference point
%     * kmax - maximal number of neighbors for each reference point
%     * Nref - number of randomly chosen reference points (n == -1 means :
%       use all points)
%     * gstart - starting value for moments
%     * gend - end value for moments
%     * past - (optional) number of samples to exclude before and after
%       each reference index, default is 0
%     * steps - (optional) number of moments to calculate, default is 32
%
%   Compute fractal dimension spectrum D(q) using moments of neighbor
%   distances for time-delay reconstructed timeseries s.
%
%   Do the main job - computing nearest neighbors for reference points.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(6,8)

if nargin < 7
	past = 0
end
if nargin < 8
	steps = 32;
end

if steps < 1
	error('Number of steps must be positive')
end

if (gend-gstart > 0) & (steps > 1)
	gammas = linspace(gstart, gend, steps);
else
	gammas = gstart;		
end

N = dlens(s,1);

try 
	atria = optparams(s, 1);
	names = fieldnames(atria);
catch
	atria = nn_prepare(data(s), 'euclidian');
	s = setoptparams(s, 1, atria);
end

% Do the main job - computing nearest neighbors for reference points 
[nn, dist] = nn_search(data(s), atria, randref(1,N,Nref), kmax, past);

out = gendimest(dist, gammas, kmin, kmin, kmax);

rs = signal(core(out(:)), s);	

a = achse(unit, 1 - (gammas(:) ./ out));
a = setname(a, 'q');
rs = setaxis(rs, 1, a);
rs = addhistory(rs,  ['Computed fractal dimension spectrum']);
rs = addcommandlines(rs, 's = fracdims(s', kmin, kmax, Nref, gstart, gend, past, steps);
rs = setyname(rs, 'D(q)');
rs = setlabel(rs, 'Renyi dimension spectrum');

end