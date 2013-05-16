%tstoolbox/mex/amutual
%   Fast, but crude auto mutual information of a scalar timeseries for the
%   timelags from zero to maxtau. The input time series should be much
%   longer than maximal timelag maxtau. The algorithm uses equidistant
%   histogram boxes, so results are bad in a mathematical sense. However,
%   a fast algorithm based on ternary search trees to store only nonempty
%   boxes is used.
%
%   Syntax:
%
%     * a = amutual(ts, maxtau, partitions)
%
%   Input arguments:
%
%     * ts - vector holding time series data
%     * maxtau - maximal time lag
%     * partitions - number of partitions for the one-dimensional
%       histogram
%
%   Output arguments:
%
%     * a - vector of length maxtau+1, holding auto mutual information
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

