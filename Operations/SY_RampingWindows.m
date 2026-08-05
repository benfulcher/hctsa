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
% ---OUTPUTS: Kendall's rank correlation tau (and its p-value) between
% segment index and each windowed statistic, for the mean, variance,
% standardized skewness/kurtosis, and AC1. Kendall's tau is used rather
% than a linear-regression slope since it is scale-invariant (comparable
% across statistics with different units) and detects any monotonic
% trend, not just a linear one -- matching "ramping" more directly than a
% slope would. Skewness/kurtosis are standardized in the usual sense
% (normalized by std^3/std^4 respectively, i.e., Matlab's skewness/
% kurtosis, not DN_Moments' convention of normalizing by std^1 regardless
% of moment order) so that a shape trend isn't conflated with the
% (separately tracked) scale trend in segVar.

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
for i = 1:numSeg
    segAC1(i) = CO_AutoCorr(z(:, i), 1, 'Fourier');
end

% ------------------------------------------------------------------------------
%% Kendall's tau (and p-value) for each statistic against segment index
% ------------------------------------------------------------------------------
segIdx = (1:numSeg)';

[out.mean_tau, out.mean_p] = SUB_kendallTrend(segIdx, segMean);
[out.var_tau, out.var_p] = SUB_kendallTrend(segIdx, segVar);
[out.skew_tau, out.skew_p] = SUB_kendallTrend(segIdx, segSkew);
[out.kurt_tau, out.kurt_p] = SUB_kendallTrend(segIdx, segKurt);
[out.ac1_tau, out.ac1_p] = SUB_kendallTrend(segIdx, segAC1);

% ------------------------------------------------------------------------------
function [tau, pval] = SUB_kendallTrend(idx, x)
    [tau, pval] = corr(idx, x(:), 'type', 'Kendall');
end
% ------------------------------------------------------------------------------

end
