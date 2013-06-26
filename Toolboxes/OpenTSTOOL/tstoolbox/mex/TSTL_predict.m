%tstoolbox/mex/predict
%   State space based prediction using nearest neighbors. The algorithms
%   computes one or more nearest neighbors to an initial state vector. The
%   images of the nearest neighbors are used to estimate to image of the
%   initial state vector. The next iteration uses the previously computed
%   image as new initial state vector .
%
%   Syntax:
%
%     * x = predict(pointset, length, k, stepsize, mode)
%
%   Input arguments:
%
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * length - number of iterations (length of prediction)
%     * k - number of nearest neighbors
%     * stepsize - prediction stepsize, usually one
%     * mode - (optional) method to estimate image of initial state vector
%          + 0 - direct prediction, no weight is applied to neighbors
%          + 1 - direct prediction, biquadratic weight is applied to
%            neighbors
%          + 2 - integrated prediction, no weight is applied to neighbors
%          + 3 - integrated prediction, biquadratic weight is applied to
%            neighbors
%
%   Output arguments:
%
%     * x - data set as double matrix, size length by D
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


