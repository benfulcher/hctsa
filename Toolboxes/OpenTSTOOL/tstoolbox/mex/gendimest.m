%tstoolbox/mex/gendimest
%   The Renyi dimension spectrum of a points set can be estimated using
%   information about the distribution of the interpoint distances. Since
%   we are interested in the scaling behaviour of the Renyi information
%   for small distances, we don't need to compute all interpoint
%   distances, the distances to k nearest neighbors for each reference
%   point are sufficient .
%
%   Robust estimation is used instead of mean square error fitting.
%
%   Syntax:
%
%     * [dimensions, moments] = gendimest(dists, gammas, kmin_low,
%       kmin_high, kmax)
%
%   Input arguments:
%
%     * dists - a matrix of size R by k which contains distances from
%       reference points to their k nearest neighbors, sorted in
%       increasing order. This matrix can be obtained by calling nn_search
%       (cf. Section ) or fnearneigh (cf. Section ) on the point set whose
%       dimension spectrum is to be investigated.
%     * gammas - vector of the moment orders
%     * kmin_low - first kmin, 1 £ kmin_low
%     * kmin_high - last kmin, kmin_low £ kmin_high < kmax
%     * kmax - highest neigbor order up to which, kmax £ k
%
%   Output arguments:
%
%     * dimensions - matrix of size length(gammas) by
%       kmin_upper-kmin_lower+1, holding the dimension estimates
%     * moments (optional) - matrix of size k by length(gammas), storing
%       the computed moments of the neigbor distances
%
%   Example:
%
%x = chaosys(25000, 0.025, [0.1 -0.1 0.02], 0);  % generate data from Lorenz sys
%tem
%x = x(5001:end,:);      % discard first 5000 samples due to transient
%[nn, dist] = fnearneigh(x, randref(1, 20000, 1000), 128, 0);
%gammas = -5:0.5:5;
%gedims = gendimest(dist, gammas, 8, 8, 128);
%plot(1-gammas./gedims', gedims)
%xlabel('q');ylabel('D_q');title('Renyi dimension')
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


