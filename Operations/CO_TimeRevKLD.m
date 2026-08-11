function out = CO_TimeRevKLD(y, tau, m, k, theilerWin, maxN)
% CO_TimeRevKLD   Kullback-Leibler divergence between forward and time-reversed embeddings.
%
% Embeds the time series in m dimensions at time delay tau (e.g., the pair
% (x_t,x_{t+tau}) for m=2, or the triple (x_t,x_{t+tau},x_{t+2tau}) for
% m=3), and estimates the Kullback-Leibler divergence between the
% distribution of these embedded points and the distribution of the same
% points with their coordinate order reversed (equivalent to embedding
% the time-reversed series). For a (statistically) time-reversible
% process, these two distributions coincide and the divergence is zero in
% the population; departures reflect time-irreversibility.
%
% cf. CO_trev/CO_TC3, which probe irreversibility via a single third-moment
% statistic of lagged pairs/triples; this is the natural full-density
% generalization (in the same sense that CO_JointNonGaussianity generalizes
% DN_Moments' skewness/kurtosis to the whole joint embedding shape), able to
% detect any asymmetry between the forward and reversed distributions, not
% just a third-moment one.
%
% cf. also Diks, van Houwelingen, Takens & DeGoede (1995), "Reversibility as
% a criterion for discriminating time series", Phys. Lett. A 201(4-5) 221,
% who compare forward and reverse embedding distributions via a different
% (U-statistic-based) route; the approach here instead estimates the
% Kullback-Leibler divergence directly using the k-NN estimator of
% Q. Wang, S.R. Kulkarni & S. Verdu, "Divergence Estimation for
% Multidimensional Densities via k-Nearest-Neighbor Distances", IEEE Trans.
% Inf. Theory 55(5) 2392 (2009), which needs no binning or KDE grid and
% remains usable at hctsa-scale sample sizes in m = 2 or 3 dimensions.
%
% The k-NN search follows the same over-fetch-then-Theiler-filter pattern as
% NL_LocalDensity: candidate neighbors are fetched via a KD-tree and any
% within a Theiler window of the query point's time index are discarded (a
% point and its temporal neighbors are strongly autocorrelated, so treating
% them as informative near-neighbors would understate the local spread);
% the rare point left with too few valid candidates falls back to an exact
% brute-force search.
%
% NOTE ON DIRECTIONALITY: KL(P||Q) and KL(Q||P) are generally different
% quantities, but *not* here: the reversed embedding Q is built by flipping
% the coordinate order of each row of P, an isometry (it preserves all
% pairwise distances, including cross-set ones), applied identically at
% every matching time index. Any purely distance-based two-sample
% divergence estimator -- like the k-NN one used here -- is therefore
% forced to return numerically identical values for KL(P||Q) and KL(Q||P)
% (verified to match to machine precision on synthetic test series); only
% one direction is computed.
%
% NOTE ON SIGNIFICANCE: only a raw divergence estimate is returned, not a
% p-value -- consecutive embedded points overlap in m-1 coordinates and are
% not independent, the same reason CO_trev/CO_TC3/CO_JointNonGaussianity
% report raw statistics only. For significance testing against a null that
% respects the series' own autocorrelation structure, compare this
% statistic to its distribution over surrogates (cf. SD_MakeSurrogates,
% SD_SurrogateTest).
%
% ---INPUTS:
% y, the input time series
%
% tau, the time delay for the embedding (can be 'ac' or 'mi', or an
%      integer, cf. BF_Embed). Default: 'ac'.
%
% m, the embedding dimension (an integer; cf. BF_Embed). Default: 2, for
%    the pairwise joint distribution (x_t,x_{t+tau}); set to 3 for the
%    triple-wise joint distribution (x_t,x_{t+tau},x_{t+2tau}).
%
% k, the number of nearest neighbors used by the k-NN divergence estimator.
%    Default: 3 (matches NL_LocalDensity's default).
%
% theilerWin, the number of temporally-adjacent points excluded from both
%             the within-set and cross-set neighbor searches (|i-j| <=
%             theilerWin), applied at matching time indices in both the
%             forward and reversed embeddings. Default: 1.
%
% maxN, the maximum number of embedded points used. The k-NN searches are
%       KD-tree-based, not the O(N^2) cost of CO_JointNonGaussianity's
%       skewness statistic (which needs this kind of cap): empirically,
%       100000 points takes ~0.3s for smooth data and ~1.6s even for
%       adversarial duplicate-heavy data that maximally triggers the
%       brute-force Theiler-window fallback (cf. local_theiler_kth). A
%       warning is issued whenever cropping actually happens. Default:
%       'full' (no cropping); set to an integer to cap runtime on
%       unusually long series.
%
% ---OUTPUTS:
% The raw k-NN estimate of KL(forward || reversed) and its magnitude. A
% unitless departure-from-reversibility measure, zero up to estimation
% noise for a reversible process (the k-NN estimator can dip slightly
% negative near a true value of zero -- expected behavior of this
% nonparametric estimator, not a bug), with no attached significance level
% (see note above).

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 'ac';
end
if nargin < 3 || isempty(m)
    m = 2;
end
if nargin < 4 || isempty(k)
    k = 3;
end
if nargin < 5 || isempty(theilerWin)
    theilerWin = 1;
end
if nargin < 6 || isempty(maxN)
    maxN = 'full';
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
Y = BF_Embed(y, tau, m, false);
if isscalar(Y) && isnan(Y) % embedding failed
    warning('Embedding failed');
    out = NaN; return
end
[Nemb, d] = size(Y);

minN = max(50, 10 * (k + theilerWin));
if Nemb < minN
    warning('Too few embedded points (%u) for a meaningful time-reversal KLD estimate at m = %u', Nemb, d);
    out = NaN; return
end

if ischar(maxN) && strcmp(maxN, 'full')
    % no cropping
elseif Nemb > maxN
    warning(['Cropping to the first %u of %u embedded points ' ...
             '(runtime cap; raise maxN or use ''full'' to disable)'], maxN, Nemb);
    Y = Y(1:maxN, :);
    Nemb = maxN;
end

% ------------------------------------------------------------------------------
%% Forward and time-reversed embeddings
% ------------------------------------------------------------------------------
% Reversing the coordinate order of each embedded row is equivalent to
% embedding the time-reversed series (up to the boundary points, which this
% avoids re-deriving from scratch):
P = Y;
Q = fliplr(Y);

% ------------------------------------------------------------------------------
%% k-NN estimate of KL(P||Q)
% ------------------------------------------------------------------------------
% (KL(Q||P) is numerically identical -- see NOTE ON DIRECTIONALITY above --
% so only one direction is computed.)
out.raw = local_kNN_KLD(P, Q, k, theilerWin);
out.abs = abs(out.raw);

end

% ------------------------------------------------------------------------------
function kld = local_kNN_KLD(A, B, k, theilerWin)
% Wang-Kulkarni-Verdu (2009) k-NN estimator of KL(P_A||P_B), where A and B
% are equal-sized (n x d) point sets whose rows are in temporal
% correspondence (row i of A and row i of B share the same underlying time
% index), so a Theiler window can be applied consistently in both searches.

n = size(A, 1);
d = size(A, 2);
kFetch = min(n - 1, k + 2 * theilerWin + 5);

% Self-search within A (for r_k(x_i), the k-th NN of x_i within A\{x_i}):
[idxA, distA] = knnsearch(A, A, 'K', kFetch + 1);
rk = local_theiler_kth(idxA, distA, k, theilerWin, A, A);

% Cross-search from A into B (for s_k(x_i), the k-th NN of x_i within B):
[idxB, distB] = knnsearch(B, A, 'K', kFetch);
sk = local_theiler_kth(idxB, distB, k, theilerWin, B, A);

good = isfinite(rk) & isfinite(sk) & rk > 0 & sk > 0;
if sum(good) < 0.5 * n
    kld = NaN; return
end
kld = (d / sum(good)) * sum(log(sk(good) ./ rk(good))) + log(n / (n - 1));

end

% ------------------------------------------------------------------------------
function kth = local_theiler_kth(idx, dist, k, theilerWin, refSet, querySet)
% For each query point i, returns the k-th nearest-neighbor distance within
% refSet, excluding any candidate with |i-j| <= theilerWin. Falls back to an
% exact brute-force search for any point whose over-fetched candidate list
% (idx, dist) didn't leave at least k valid (non-excluded) neighbors.

n = size(idx, 1);
nRef = size(refSet, 1);
kth = nan(n, 1);
for i = 1:n
    validDists = dist(i, abs(idx(i, :) - i) > theilerWin);
    if length(validDists) >= k
        kth(i) = validDists(k);
    else
        allDists = sqrt(sum((refSet - querySet(i, :)).^2, 2));
        allDists(abs((1:nRef)' - i) <= theilerWin) = Inf;
        sortedDists = sort(allDists);
        if sortedDists(k) < Inf
            kth(i) = sortedDists(k);
        end
    end
end

end
