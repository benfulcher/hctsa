%tstoolbox/mex/cao
%   This mex-file applies Cao's method to the input data set. If the data
%   set contains points of dimension D, it computes E and E* for a data
%   set of dimension 1 (taken from the first column of the input data
%   set), then for a data set of dimension 2 (taken from the first two
%   columns) up to a dimension of D. Optionally, this algorithm extends
%   Cao's method in a straightforward manner to use more than one nearest
%   neighbors.
%
%   Syntax:
%
%     * [E, E*] = cao(pointset, query_indices, k)
%
%   Input arguments:
%
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * query_indices - query points are taken out of the pointset,
%       query_indices is a vector of length R which contains the indices
%       of the query points (indices may vary from 1 to N)
%     * k - number of nearest neighbors to compute. Cao's method can be
%       extended to use more than only the first nearest neighbor (k=1).
%
%   Output arguments:
%
%     * E and E* are vectors of size D. Please refer the Cao's article for
%       a precise description of their meaning.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

