function out = SY_RampingWindows(y, numSeg)
% SY_RampingWindows    Monotonic trend ('ramping') in windowed statistics
%
% Splits the time series into numSeg non-overlapping segments, computes
% the mean, variance, skewness, kurtosis, and lag-1 autocorrelation (AC1)
% within each segment, and quantifies whether each of these quantities
% trends monotonically across the segments (e.g., a variance that ramps
% up steadily across the series, rather than merely fluctuating).
%
% Existing hctsa stationarity operations (SY_SlidingWindow, SY_DriftingMean,
% SY_StatAv, SY_LocalDistributions) summarize the *spread* of windowed
% statistics (std, range, entropy) -- order-agnostic measures that don't
% distinguish a monotonic ramp from random fluctuation of the same
% magnitude. This operation targets that gap directly.
%
% ---INPUTS:
% y, the input time series
%
% numSeg, the number of non-overlapping segments to divide the time series
%         into (default: 10). Non-overlapping segments are used
%         deliberately (rather than SY_SlidingWindow's overlapping
%         windows): overlap between windows would induce artificial
%         serial correlation between adjacent window-statistics, which
%         would inflate the apparent monotonic trend independent of any
%         real ramping in the data.
%
% ---OUTPUTS: Kendall's rank correlation tau, and Pearson's linear
% correlation r (each with its p-value), between segment index and each
% windowed statistic, for the mean, variance, standardized
% skewness/kurtosis, AC1, and asymAC1 (see below). Kendall's tau is
% scale-invariant and detects any monotonic trend, not just a linear one
% -- matching "ramping" more directly than a slope would, and more
% robust to outlying segments. Pearson's r is included alongside it as a
% complementary, effect-size-like measure of specifically *linear*
% ramping (r^2 is the fraction of across-segment variance explained by a
% linear trend); expect the two to agree closely for a clean linear ramp
% and diverge for a monotonic-but-nonlinear one (e.g. a ramp that
% plateaus). Skewness/kurtosis are standardized in the usual sense
% (normalized by std^3/std^4 respectively, i.e., Matlab's skewness/
% kurtosis, not DN_Moments' convention of normalizing by std^1 regardless
% of moment order) so that a shape trend isn't conflated with the
% (separately tracked) scale trend in segVar.
%
% asymAC1 is a nonlinear, asymmetric variant of the lag-1 autocorrelation:
% mean(x_t * x_{t+1}^2), with x z-scored *within* each segment (unlike the
% other statistics above, this one needs the explicit per-segment
% z-scoring, since it's neither scale-invariant like AC1/skewness/
% kurtosis nor a raw moment tracked deliberately like mean/variance).
% It's zero in expectation for any time-reversible linear process (same
% spirit as CO_trev's third-moment reversibility statistic, computed
% densely at a single lag over the whole series rather than trended
% across segments), so a ramp in asymAC1 flags a trend specifically in
% the series' local time-asymmetry/nonlinearity, distinct from a trend in
% any of the other (symmetric) segment statistics above.

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
N = length(y);

if nargin < 2 || isempty(numSeg)
    numSeg = 10;
end

minNumSeg = 5; % need enough segments for a meaningful trend statistic
minSegLength = 20; % heuristic minimum for meaningful skewness/kurtosis/AC1 estimates

if numSeg < minNumSeg
    error('numSeg = %u is too few segments for a meaningful trend statistic (need >= %u)', numSeg, minNumSeg);
end

segLength = floor(N / numSeg);
if segLength < minSegLength
    warning('Time series (N = %u) too short for %u segments of a meaningful length', N, numSeg);
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Segment the time series (non-overlapping, discarding any remainder)
% ------------------------------------------------------------------------------
z = reshape(y(1:segLength * numSeg), segLength, numSeg); % segLength x numSeg

% ------------------------------------------------------------------------------
%% Within-segment statistics
% ------------------------------------------------------------------------------
segMean = mean(z);
segVar = var(z);
% Standardized (not DN_Moments-style raw/std) skewness and kurtosis: keeps
% the shape-trend signal from being conflated with the (separately
% tracked) scale-trend signal in segVar, since DN_Moments normalizes by
% std^1 regardless of moment order rather than std^3/std^4.
segSkew = skewness(z);
segKurt = kurtosis(z);
segAC1 = zeros(1, numSeg);
segAsymAC1 = zeros(1, numSeg);
for i = 1:numSeg
    segAC1(i) = CO_AutoCorr(z(:, i), 1, 'Fourier');
    zseg = zscore(z(:, i)); % z-scored *within* this segment
    segAsymAC1(i) = mean(zseg(1:end - 1) .* zseg(2:end).^2);
end

% ------------------------------------------------------------------------------
%% Kendall's tau and Pearson's r (each with p-value) against segment index
% ------------------------------------------------------------------------------
segIdx = (1:numSeg)';

[out.mean_tau, out.mean_p] = SUB_trend(segIdx, segMean, 'Kendall');
[out.var_tau, out.var_p] = SUB_trend(segIdx, segVar, 'Kendall');
[out.skew_tau, out.skew_p] = SUB_trend(segIdx, segSkew, 'Kendall');
[out.kurt_tau, out.kurt_p] = SUB_trend(segIdx, segKurt, 'Kendall');
[out.ac1_tau, out.ac1_p] = SUB_trend(segIdx, segAC1, 'Kendall');
[out.asymac1_tau, out.asymac1_p] = SUB_trend(segIdx, segAsymAC1, 'Kendall');

[out.mean_pearson_r, out.mean_pearson_p] = SUB_trend(segIdx, segMean, 'Pearson');
[out.var_pearson_r, out.var_pearson_p] = SUB_trend(segIdx, segVar, 'Pearson');
[out.skew_pearson_r, out.skew_pearson_p] = SUB_trend(segIdx, segSkew, 'Pearson');
[out.kurt_pearson_r, out.kurt_pearson_p] = SUB_trend(segIdx, segKurt, 'Pearson');
[out.ac1_pearson_r, out.ac1_pearson_p] = SUB_trend(segIdx, segAC1, 'Pearson');
[out.asymac1_pearson_r, out.asymac1_pearson_p] = SUB_trend(segIdx, segAsymAC1, 'Pearson');

% ------------------------------------------------------------------------------
function [rho, pval] = SUB_trend(idx, x, corrType)
    [rho, pval] = corr(idx, x(:), 'type', corrType);
end
% ------------------------------------------------------------------------------

end
