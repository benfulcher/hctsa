%tstoolbox/mex/boxcount
%   boxcount is a fast algorithm that partitions a data set of points into
%   equally spaced and sized boxes. The algorithm is based on Robert
%   Sedgewick's Ternary Search Trees which offer a fast and efficient way
%   to create and search a multidimensional histogram. Empty boxes require
%   no storage space, therefore the maximum number of boxes (and memory)
%   used can not exceed the number of points in the data set, regardless
%   of the data set's dimension and the number of partitions per axis.
%
%   During processing, data values are scaled to be within the range
%   [0,1]. All columns of the input matrix are scaled by the same factor,
%   so no skewing is introduced into the point set.
%
%   Syntax:
%
%     * [a,b,c] = boxcount(point_set, partitions)
%
%   Input arguments:
%
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D. D is limited
%       to 128.
%     * partitions - number of partitions per axis, limited to 16384. For
%       convenience, if a vector is given, boxcount will iterate over all
%       values of this vector.
%
%   Output arguments:
%
%     * a - vector of size D with: log2(sum(Number of nonempty boxes))
%     * b - vector of size D with: sum(p * log2(p)) , where p is the
%       relative frequency of points falling into a box
%     * c - vector of size D with: log2(sum(p*p)), where p is the relative
%       frequency of points falling into a box
%
%   Example:
%
%p = rand(50000, 4);
%p = p - min(min(p));
%p = p ./ max(max(p));
%[a,b,c] = boxcount(p, 16)
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

