%tstoolbox/mex/range_search
%   Syntax:
%
%     * [count, neighbors] = range_search(pointset, atria, query_points,
%       r)
%     * [count, neighbors] = range_search(pointset, atria, query_indices,
%       r, exclude)
%
%   Input arguments:
%
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * atria - output of (cf. Section )nn_prepare for pointset
%     * query_points - a R by D double matrix containing the coordinates
%       of the query points, organized as R points of dimension D
%     * query_indices - query points are taken out of the pointset,
%       query_indices is a vector of length R which contains the indices
%       of the query points
%     * r - range or search radius (r > 0)
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. For example if the index of the query point is 124 and
%       exclude is set to 3, points with indices 121 to 127 are omitted
%       from search. Using exclude = 0 means: exclude self-matches
%
%   Output arguments:
%
%     * count - a vector of length R contains the number of points within
%       distance r to the corresponding query point
%     * neighbors - a Matlab cell structure of size R by 2 which contains
%       vectors of indices and vectors of distances to the neighbors for
%       each given query point. This output argument can not be stored in
%       a standard Matlab matrix because the number of neighbors within
%       distance r is not the same for all query points. The vectors if
%       indices and distances for one query point have exactly the length
%       that is given in count. The values in the distances vectors are
%       not sorted..
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


