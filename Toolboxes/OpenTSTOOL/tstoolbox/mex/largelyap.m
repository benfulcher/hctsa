%tstoolbox/mex/largelyap
%   largelyap is an algorithm very similar to the Wolf algorithm , it
%   computes the average exponential growth of the distance of neighboring
%   orbits via the prediction error. The increase of the prediction error
%   vs the prediction time allows an estimation of the largest lyapunov
%   exponent.
%
%   Syntax:
%
%     * x = largelyap(pointset, query_indices, taumax, k exclude)
%     * x = largelyap(atria, pointset, query_indices, taumax, k exclude)
%
%   Input arguments:
%
%     * atria - output of nn_prepare for pointset (optional)
%       (cf. Section )
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * query_indices - query points are taken out of the pointset,
%       query_indices is a vector of length R which contains the indices
%       of the query points (indices may vary from 1 to N)
%     * taumax - maximal time shift
%     * k - number of nearest neighbors to compute
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. For example if the index of the query point is 124 and
%       exclude is set to 3, points with indices 121 to 127 are omitted
%       from search. Using exclude = 0 means: exclude self-matches
%
%   Output arguments:
%
%     * x - vector of length taumax+1, x(tau) = 1/Nref *
%       sum(log2(dist(reference point + tau, nearest neighbor +
%       tau)/dist(reference point, nearest neighbor)))
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


