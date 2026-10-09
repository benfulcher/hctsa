function out = SC_HiguchiFD(y, kmax)
% SC_HiguchiFD   Higuchi's fractal dimension of the waveform of a time series.
%
% Treats the series as a curve in the (index, value) plane and measures its
% length at a range of scales by Higuchi's curve-length method. At each scale k
% the series is split into k interleaved subsequences (every k-th point,
% starting from each of the first k points); the length of each subsequence's
% curve is measured, normalized for its number of points, and the k lengths are
% averaged to give L(k). For a fractal curve L(k) ~ k^(-D), so the fractal
% dimension D is the slope of log L(k) against log(1/k), found by a robust
% linear fit. Smooth series give D near 1, Brownian motion about 1.5, and white
% noise near 2.
%
% ---INPUTS:
% y, the input time series
% kmax, the maximum interleaving scale to include in the fit (default:
%       min(8, max(5, floor(N/10))), a small, fixed value; see NOTES for why a
%       large kmax is harmful)
%
% ---OUTPUTS:
% HFD, the Higuchi fractal dimension estimate (the log-log slope)
% intercept, the intercept of the fit
% se_HFD, the standard error of the slope
% ssr, the mean squared residual of the fit
% resac1, the lag-1 autocorrelation of the residuals (large values indicate
%       systematic curvature in the log-log plot, so that a single power law
%       describes the scaling poorly)
% (A scalar NaN is returned if fewer than 5 scales give a usable curve length.)
%
% ---REFERENCES:
% T. Higuchi, "Approach to an irregular time series on the basis of the fractal
% theory", Physica D 31(2), 277-283 (1988).
%
% ---NOTES:
% Unlike DFA-style estimators (e.g., SC_FluctAnal), which detrend a
% cumulative-sum (integrated) profile within windows and measure how fluctuations
% scale with window size (targeting long-range-correlation-type self-affinity,
% e.g., the Hurst exponent of fractional Brownian motion), Higuchi's method
% measures the curve length of the raw series directly, across interleaved
% subsequences, with no integration or detrending step. It is more sensitive to
% short-range waveform roughness than to long-range correlation structure. It is
% also distinct from hctsa's embedding-based dimension estimators
% (NL_FractalDimensions, NL_Dimensions), which measure the geometry of a
% time-delay-embedded attractor rather than the waveform of the raw series.
%
% kmax must be kept small and not scaled up with series length N. Tested against
% known ground truth: for white noise (D = 2) and Brownian motion (D = 1.5) the
% fitted dimension is robust to kmax across 8-200. But for a smooth periodic
% signal (D about 1) the fit is essentially exact at kmax = 8 (D = 1.011) and
% degrades badly as kmax grows relative to the signal's characteristic timescale
% (D = 1.76 at kmax = 200, on a period-137 sinusoid): once the interleaving scale
% exceeds roughly half the oscillation period, successive sampled points alias
% across cycles, inflating the apparent curve length at large k and biasing the
% fitted slope upward. The default matches common practice in the Higuchi-FD
% literature (a small, fixed kmax, typically about 8-16 regardless of N).

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
if size(y, 2) > size(y, 1)
	y = y'; % ensure a column vector
end
N = length(y);

if nargin < 2 || isempty(kmax)
	kmax = min(8, max(5, floor(N / 10))); % small and fixed once N is not tiny -- see NOTES
end

% ------------------------------------------------------------------------------
%% Compute the curve length L(k) at each scale k = 1,...,kmax
% ------------------------------------------------------------------------------
logL = nan(kmax, 1);
logInvk = nan(kmax, 1);
for k = 1:kmax
	Lk = nan(k, 1);
	for m = 1:k
		nMax = floor((N - m) / k);
		if nMax < 1
			continue
		end
		idx = m:k:(m + nMax * k);
		Lk(m) = sum(abs(diff(y(idx)))) * (N - 1) / (nMax * k) / k;
	end
	Lbar = mean(Lk,'omitnan');
	if ~isnan(Lbar) && Lbar > 0
		logL(k) = log(Lbar);
		logInvk(k) = log(1 / k);
	end
end

goodK = ~isnan(logL);
if sum(goodK) < 5
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Robust linear fit of log(L(k)) against log(1/k)
% ------------------------------------------------------------------------------
[linfit, stats] = robustfit(logInvk(goodK), logL(goodK));

out.HFD = linfit(2); % the Higuchi fractal dimension estimate (log-log slope)
out.intercept = linfit(1);
out.se_HFD = stats.se(2); % standard error on the dimension estimate
out.ssr = mean(stats.resid.^2); % mean squared residual of the linear fit
out.resac1 = CO_AutoCorr(stats.resid, 1, 'Fourier'); % residual autocorrelation
% (large resac1 indicates systematic curvature in the log-log plot -- a
% single power law is a poor description of the scaling behavior)

end
