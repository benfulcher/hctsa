%tstoolbox/mex/corrsum
%   The topics correlation sum and correlation dimension estimation can
%   also be found here.
%
%   Syntax:
%
%     * [c, d] = corrsum(pointset, query_indices, range, exclude)
%     * [c, d] = corrsum(pointset, query_indices, range, exclude, bins)
%     * [c, d] = corrsum(atria, pointset, query_indices, range, exclude)
%     * [c, d] = corrsum(atria, pointset, query_indices, range, exclude,
%       bins)
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
%     * range - search range, may be given in one of two ways
%          + If only a single value is given, this value is taken as
%            maximal search radius relative to the attractor diameter (0 <
%            relative_range < 1). The minimal search radius is determined
%            automatically be searching for the minimal interpoint
%            distance in the data set.
%          + If a vector of length two is given, the values are
%            interpreted as absolut minimal and maximal search radius.
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. For example if the index of the query point is 124 and
%       exclude is set to 3, points with indices 121 to 127 are omitted
%       from search. Using exclude = 0 means: exclude self-matches
%     * bins - number of distance values at which the correlation sum is
%       evaluated, defaults to 32
%
%   Output arguments:
%
%     * c - vector of correlation sums, length(c) = bins
%     * d - vector of the corresponding distances at which the correlation
%       sums (stored in c) were computed. d is exponentially spaced,
%       length(c) = bins
%
%   Example:
%
%x = chaosys(25000, 0.025, [0.1 -0.1 0.02], 0);  % generate data from Lorenz sys
%tem
%x = x(5001:end,:);      % discard first 5000 samples due to transient
%% now compute correlation sum up to five percent of attractor diameter
%[c,d] = corrsum(x, randref(1,20000, 1000), 0.05, 0);
%loglog(d,c)     % and show the result as log-log plot
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

