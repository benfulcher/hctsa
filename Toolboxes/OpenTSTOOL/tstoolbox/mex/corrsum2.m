%tstoolbox/mex/corrsum2
%   This is an extended version of the correlation sum algorithm. It tries
%   to accelerate the computation of the correlation sum by using a
%   different number of reference points at each length scale. For large
%   length scales, only a few number of reference points will be used
%   since for this scale, quite a lot of neighbors will fall within this
%   range (and also the search time will be high). The smaller the length
%   scale, the more reference points are used. The algorithm tries to keep
%   the number of pairs found within each range roughly constant at Npairs
%   to ensure a good statistic even for the smallest length scales.
%   However, the number of reference points actually used may be limited
%   to be within [Nref_min Nref_max] to give at least some control to the
%   user. All reference points are chosen randomly from the data set
%   without reoccurences of the same index.
%
%   Syntax:
%
%     * [c, d, e, f, g] = corrsum(pointset, Npairs, range, exclude)
%     * [c, d, e, f, g] = corrsum(pointset, Npairs, range, exclude, bins)
%     * [c, d, e, f, g] = corrsum(pointset, Npairs, range, exclude, bins,
%       opt_flag)
%     * [c, d, e, f, g] = corrsum(atria, pointset, Npairs, range, exclude)
%     * [c, d, e, f, g] = corrsum(atria, pointset, Npairs, range, exclude,
%       bins)
%     * [c, d, e, f, g] = corrsum(atria, pointset, Npairs, range, exclude,
%       bins, opt_flag)
%
%   Input arguments:
%
%     * atria - output of nn_prepare for pointset (optional)
%       (cf. Section )
%     * pointset - a N by D double matrix containing the coordinates of
%       the point set, organized as N points of dimension D
%     * Npairs - Number of pairs to find within each length scale. The
%       algorithm will adapt the number of reference points while
%       computing the correlation sum. Reference points are chosen
%       randomly from the pointset. Optionally, a vector of the form
%       [Npairs Nref_min Nref_max] may be given. For no length scale less
%       than Nref_min reference points will be used. Additionally, not
%       more than Nref_max reference points will be used at all.
%     * range - search range, may be given in one of two ways
%          + If only a single value is given, this value is taken as
%            maximal search radius relative to attractor diameter (0 <
%            relative_range < 1). The minimal search radius is determined
%            automatically be searching for the minimal interpoint
%            distance in the data set.
%          + If a vector of length two is given, the values are
%            interpreted as absolut minimal and maximal search radius.
%     * exclude - in case the query points are taken out of the pointset,
%       exclude specifies a range of indices which are omitted from
%       search. E.g. if the index of the query point is 124 and exclude is
%       set to 3, points with indices 121 to 127 are omitted from search.
%       exclude = 0 means : exclude self-matches
%     * bins - number of distance values at which the correlation sum is
%       evaluated, defaults to 32
%     * opt_flag - optional flag to control the algorithm:
%          + 0 - Use euclidian distance, be verbose, don't allow to count
%            a pair of points twice
%          + 1 - Use maximum distance, be verbose, don't allow to count a
%            pair of points twice
%          + 2 - Use euclidian distance, be verbose, allow to count a pair
%            of points twice
%          + 3 - Use maximum distance, be verbose, allow to count a pair
%            of points twice
%          + 4 - Use euclidian distance, be silent, don't allow to count a
%            pair of points twice
%          + 5 - Use maximum distance, be silent, don't allow to count a
%            pair of points twice
%          + 6 - Use euclidian distance, be silent, allow to count a pair
%            of points twice
%          + 7 - Use maximum distance, be silent, allow to count a pair of
%            points twice
%       If the preprocessing output atria is given, the type of metric
%       used to create this overrides the settings by opt_flag.
%
%   Output arguments:
%
%     * c - vector of correlation sums, length(c) = bins
%     * d - vector of the corresponding distances at which the correlation
%       sums (stored in c) where computed. d is exponentially spaced,
%       length(c) = bins
%     * e - vector of the number of pairs found within this range,
%       length(e) = bins
%     * f - vector of the number of total pairs that were tested,
%       length(f) = bins
%     * g - vector containing the indices of the reference points actually
%       used by the algorithm.
%
%   Example:
%
%x = chaosys(25000, 0.025,  [0.1 -0.1 0.02], 0);
%x = x(5001:end,:);              % discard first 5000 samples due to transient
%% now compute correlation sum up to five percent of attractor diameter
%[c,d] = corrsum2(x,[1000 100 2000], 0.05, 200);
%loglog(d,c)        % and show the result as log-log plot
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


