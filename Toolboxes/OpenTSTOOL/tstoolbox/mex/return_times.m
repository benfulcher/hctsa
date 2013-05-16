%tstoolbox/mex/return_times
%   return_time may be used to find hidden periodicity in multivariate
%   data, e.g. embedded time series data. It computes a histogram of
%   return times. For any given reference point, return_time calculates
%   the time span until the time series returns to that location in phase
%   space (by means of nearest neighbors). A histogram of these time spans
%   is computed. Strong peaks in this histogram might be a sign of
%   periodicity in the data.
%
%   Syntax:
%
%     * r = return_time(pointset, query_indices, k, max_time, exclude)
%     * r = return_time(atria, pointset, query_indices, k, max_time,
%       exclude)
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
%     * k - number of nearest neighbors to be determined
%     * max_time - integer scalar, gives an upper limit for return times
%       that should be considered.
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. For example if the index of the query point is 124 and
%       exclude is set to 3, points with indices 121 to 127 are omitted
%       from search. Using exclude = 0 means: exclude self-matches
%
%   Output arguments:
%
%     * r - vector of length max_time, containing the histogram of return
%       times
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


