%tstoolbox/mex/nn_prepare
%   The intention of this mex-file was to reduce the computational
%   overhead of preprocessing for nearest neighbor or range searching.
%   With nn_prepare it is possible to do the preprocessing for a given
%   point set only once and save the created tree structure into a Matlab
%   variable. This Matlab variable, usually called atria, can then be used
%   for repeated neighbor searches on the same point set. Most mex-files
%   that rely on nearest neighbor or range search offer the possibility to
%   use this variable atria as optional input argument. However, if the
%   underlying point set is altered in any way, the proprocessing has to
%   be repeated for the new point set. If the preprocessing output does
%   not belong to the given point set, wrong results or program
%   termination may occur.
%
%   Syntax:
%
%     * atria = nn_prepare(pointset)
%     * atria = nn_prepare(pointset, metric)
%     * atria = nn_prepare(pointset, metric, clustersize)
%
%   Input arguments:
%
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * metric - (optional) either 'euclidian' or 'maximum' (default is
%       'euclidian')
%     * clustersize - (optional) threshold for clustering algorithm,
%       defaults to 64
%
%   Example:
%
%pointset = rand(40000, 3);
%atria = nn_prepare(pointset);
%[c, d] = corrsum(atria, pointset, 1:17:40000, 0.05, 0);
%plot(log(d), log(c))
%D = takens_estimator(atria, pointset, 1:17:40000, 0.05, 0)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


