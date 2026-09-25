function out = CO_JointNonGaussianity(y, tau, m, theilerWin, maxN)
% CO_JointNonGaussianity   Tests for non-Gaussianity of the joint, time-lagged embedding distribution.
%
% Embeds the time series in m dimensions at time delay tau (e.g., the pair
% (x_t,x_{t+tau}) for m=2, or the triple (x_t,x_{t+tau},x_{t+2tau}) for
% m=3) and tests whether the resulting point cloud is consistent with a
% multivariate Gaussian.
%
% cf. existing Gaussianity tests acting on marginal distribution
% (HT_DistributionTest, DN_CompareKSFit), which all
% A linear (e.g., AR(1)) Gaussian process has a Gaussian marginal
% *and* a Gaussian joint embedding distribution;
% a nonlinear or non-reversible process can look
% Gaussian marginally while its lagged joint distribution is visibly
% non-elliptical (curved, multimodal, or heavy/light-tailed along
% directions the marginal alone cannot see).
% cf. also time-irreversibility metrics like CO_trev/CO_TC3
% (which use a single third-moment statistic of lagged
% pairs/triples as a nonlinearity probe) but this is perhaps more general as
% : it tests the whole joint shape rather than one moment combination.
%
% Two complementary statistics are computed, both based on Mardia's
% (1970) classical multivariate normality measures, chosen because they
% generalize to any embedding dimension m via the same formula (so m=2
% and m=3 are the same code path) and because they reduce, at m=1, to
% ordinary skewness/kurtosis -- the natural multivariate extension of
% what DN_Moments already computes:
%   (i)  Mardia's multivariate skewness, b1 -- detects asymmetry/curvature
%        of the joint distribution (e.g., a banana-shaped point cloud).
%        Population value is 0 for any joint Gaussian.
%   (ii) Mardia's multivariate kurtosis, b2 -- detects joint tail weight/
%        peakedness relative to a Gaussian ellipsoid. Population value is
%        m(m+2) for any joint Gaussian (e.g., 8 at m=2, 15 at m=3).
% As a complementary, distribution-shape-sensitive check, the squared
% Mahalanobis distances of each embedded point to the sample mean (which
% are exactly the per-point terms underlying Mardia's kurtosis) are
% compared against their theoretical shape under joint Gaussianity,
% chi^2_m, via a Kolmogorov-Smirnov D-statistic -- this can catch
% departures (e.g., a bimodal or ring-shaped cloud) that the two summary
% moments can miss.
%
% NOTE ON SIGNIFICANCE: only *raw* statistics are returned, not p-values.
% Mardia's classical asymptotic null distributions assume the N embedded
% points are iid draws, but consecutive embedded vectors overlap in m-1
% coordinates and are therefore strongly autocorrelated, which inflates
% the naive asymptotic test statistics (empirically, up to ~30% false
% positives at a nominal 5% level on a purely linear-Gaussian AR(1)
% process, worse at higher m). This is the same reason CO_trev/CO_TC3
% report raw statistics rather than p-values; for significance
% testing against a null that respects the series' own autocorrelation
% structure, compare these statistics to their distribution over
% surrogates (cf. SD_SurrogateTest, SD_MakeSurrogates).
%
% The skewness statistic additionally excludes near-diagonal pairs
% (|i-j| <= theilerWin) from its double sum: for correlated (not just
% independent) jointly-Gaussian points, the third moment of their
% Mahalanobis inner product is *not* zero (only the independent case has
% this symmetry), so nearby, strongly-autocorrelated pairs bias the raw
% statistic away from zero even under true joint Gaussianity. (The
% same rationale as NL_RQA's Theiler window, but applied to a third-moment sum
% instead of a distance threshold).
% Empirically this removes most, but not
% all, of the bias (e.g., at m=7 on an AR(1) process, ~0.18 -> ~0.07); the
% residual is the classical small-sample bias of using the *sample*
% covariance to whiten the same points being tested (present even for iid
% data), which widening the window further does not touch.
%
% ---INPUTS:
% y, the input time series
%
% tau, the time delay for the embedding (can be 'ac' or 'mi', or an
%      integer, cf. BF_Embed). Default: 'ac'.
%
% m, the embedding dimension (can be an integer, or {'fnn',th}, cf.
%    BF_Embed). Default: 2, for the pairwise joint distribution
%    (x_t,x_{t+tau}); set to 3 for the triple-wise joint distribution
%    (x_t,x_{t+tau},x_{t+2tau}).
%
% theilerWin, the number of temporally-adjacent embedded points excluded
%             from the skewness double sum (|i-j| <= theilerWin), to
%             reduce the correlated-pair bias described above: {'ac', k}
%             for k times the first zero-crossing of the autocorrelation
%             function, or a number of samples (see BF_TheilerWindow).
%             Default: {'ac', 1}.
%
% maxN, the maximum number of embedded points used for the skewness
%       statistic (default: 10000; 'full' to disable). Legacy cap: the
%       skewness double sum used to be evaluated through an N x N Gram
%       matrix, O(N^2) in time and memory, which this bounded. It is now
%       evaluated exactly through the third-moment tensor in O(N d^3)
%       (see the code), so the cap costs nothing to lift -- it is honored
%       only so that values are unchanged from earlier computations on
%       series with more than maxN embedded points. The mean, covariance,
%       kurtosis and KS statistic always use the full embedded series.
%
% ---OUTPUTS:
% Mardia's raw multivariate skewness (Theiler-windowed) and multivariate
% kurtosis, and the Mahalanobis-distance-vs-chi^2 Kolmogorov-Smirnov
% D-statistic. All are unitless departure-from-joint-Gaussianity
% magnitudes with no attached significance level (see note above).
%
% cf. K.V. Mardia, "Measures of multivariate skewness and kurtosis with
% applications", Biometrika 57(3) 519 (1970).

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
if nargin < 4 || isempty(theilerWin)
    theilerWin = {'ac', 1};
end
theilerWin = BF_TheilerWindow(y, theilerWin);
if isnan(theilerWin) % the autocorrelation function never crosses zero
    warning('No autocorrelation zero-crossing to set the Theiler window')
    out = NaN; return
end
if nargin < 5 || isempty(maxN)
    maxN = 10000; % caps the O(N^2) skewness Gram matrix at ~800MB
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

% Need enough points to reliably estimate a d x d covariance matrix and
% for the higher-moment statistics below to be reasonably stable:
minN = max(30, 10 * d * (d + 2));
if Nemb < minN
    warning('Too few embedded points (%u) for a meaningful joint-Gaussianity test at m = %u', Nemb, d);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Center and whiten
% ------------------------------------------------------------------------------
mu = mean(Y, 1);
Yc = Y - mu;
S = (Yc' * Yc) / Nemb; % Mardia's convention: divide by N, not N-1

[L, p] = chol(S, 'lower');
if p ~= 0 || rcond(S) < 1e-10
    % Near-singular covariance -- typically means tau is too small relative
    % to the series' correlation length, so consecutive embedded
    % coordinates are nearly collinear
    warning('Embedded covariance matrix is near-singular (tau too small?)');
    out = NaN; return
end

X = L \ Yc'; % d x Nemb; column i is the whitened point L^{-1}*(Y(i,:)-mu)'
D2 = (sum(X.^2, 1))'; % squared Mahalanobis distances, Nemb x 1

% ------------------------------------------------------------------------------
%% Mardia's multivariate kurtosis (O(N), uses all embedded points)
% ------------------------------------------------------------------------------
% Raw statistic; population value under joint Gaussianity is d*(d+2)
% regardless of dimension or sample dependence structure.
out.mardiaKurt = mean(D2.^2);

% ------------------------------------------------------------------------------
%% Mahalanobis-distance-vs-chi^2_d Kolmogorov-Smirnov statistic (O(N))
% ------------------------------------------------------------------------------
x1 = sort(unique(round(D2 * 1e6) / 1e6));
if x1(1) > min(D2), x1 = [min(D2); x1]; end
if x1(end) < max(D2), x1 = [x1; max(D2)]; end
mycdf = [x1, chi2cdf(x1, d)];
[~, ~, out.mahalKSstat] = kstest(D2, 'CDF', mycdf);

% ------------------------------------------------------------------------------
%% Mardia's multivariate skewness
% ------------------------------------------------------------------------------
% b_{1,d} = mean over pairs (i,j) of (x_i' x_j)^3 for whitened x. Rather than
% forming the N x N Gram matrix (O(N^2) time and memory, which is what the
% maxN cap was guarding), use the identity
%   sum_{i,j} (x_i' x_j)^3 = sum_{a,b,c} ( sum_i x_ia x_ib x_ic )^2 = ||M||^2,
% where M is the d x d x d third-moment tensor -- O(N d^3) -- and then
% subtract the O(N w) contribution of the pairs inside the Theiler band
% |i - j| <= w, which the mean excludes. Identical result, no N^2 term.
if ~(ischar(maxN) && strcmp(maxN, 'full')) && Nemb > maxN
    % (kept for backwards compatibility of values with earlier computations;
    % no longer needed for cost)
    Xskew = X(:, 1:maxN);
else
    Xskew = X;
end
Nskew = size(Xskew, 2);
numOffBand = Nskew^2 - ((2 * theilerWin + 1) * Nskew - theilerWin * (theilerWin + 1));
if theilerWin >= Nskew || numOffBand <= 0
    warning('theilerWin too large relative to the (possibly subsampled) skewness sample size');
    out.mardiaSkew = NaN;
else
    % Full double sum via the third-moment tensor:
    Xt = Xskew'; % Nskew x d
    fullSum = 0;
    for a = 1:d
        for b = 1:d
            Mab = (Xt(:, a) .* Xt(:, b))' * Xt; % 1 x d: sum_i x_ia x_ib x_ic over c
            fullSum = fullSum + sum(Mab.^2);
        end
    end
    % Pairs inside the Theiler band (including i == j), to exclude:
    bandSum = sum(sum(Xskew.^2, 1).^3); % i == j: (x_i' x_i)^3
    for k = 1:min(theilerWin, Nskew - 1)
        ip = sum(Xskew(:, 1:end - k) .* Xskew(:, 1 + k:end), 1); % x_i' x_{i+k}
        bandSum = bandSum + 2 * sum(ip.^3);
    end
    out.mardiaSkew = (fullSum - bandSum) / numOffBand;
end

end
