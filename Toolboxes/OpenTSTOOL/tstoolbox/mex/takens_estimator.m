%tstoolbox/mex/takens_estimator
%   Syntax:
%
%     * D = takens_estimator(pointset, query_indices, relative_range,
%       exclude)
%     * D = takens_estimator(atria, pointset, query_indices,
%       relative_range, exclude)
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
%     * relative_range - search radius, relative to attractor diameter (0
%       < relative_range < 1)
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. For examples if the index of the query point is 124 and
%       exclude is set to 3, points with indices 121 to 127 are omitted
%       from search. Using exclude = 0 means: exclude self-matches
%
%   Output arguments:
%
%     * D - scalar value, estimation of correlation dimension
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


