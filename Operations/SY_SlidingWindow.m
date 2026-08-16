function out = SY_SlidingWindow(y, windowStat, acrossWinStat, numSeg, incMove)
% SY_SlidingWindow  Sliding window measures of stationarity.
%
% This function is based on sliding a window along the time series, measuring
% some quantity in each window, and outputting some summary of this set of local
% estimates of that quantity.
%
% Another way of saying it: calculate 'windowStat' in each window, and computes
% 'acrossWinStat' for the set of statistics calculated in each window.
%
% ---INPUTS:
%
% y, the input time series
%
% windowStat, the measure to calculate in each window:
%               (i) 'mean', mean
%               (ii) 'std', standard deviation
%               (iii) 'ent', distribution entropy
%               (iv) 'mom3', skewness
%               (v) 'mom4', kurtosis
%               (vi) 'mom5', the fifth moment of the distribution
%               (vii) 'lillie', the Lilliefors test statistic (KS-type distance
%                   from a Gaussian CDF) -- not its p-value, which MATLAB clips
%                   to a lookup-table boundary once data is clearly non-Gaussian
%               (viii) 'AC1', the lag-1 autocorrelation
%               (ix) 'permen', normalized Permutation Entropy, PermEn(3,1)
%                   -- a single ordinal/dynamical-complexity entropy measure,
%                   used in place of ApEn/SampEn: those are biased on the
%                   short windows this function typically operates on
%                   (numSeg=10 segments), whereas permutation entropy's
%                   pattern-count estimator degrades more gracefully at
%                   small sample sizes
%               (x) 'specen', normalized Shannon spectral entropy of the
%                   window's power spectrum (window mean removed before the
%                   FFT) -- captures whether the *frequency content* drifts
%                   across the record, distinct from all the above
%                   time-domain measures
%               (xi) 'asymAC1', mean(x_t*x_{t+1}^2) with x z-scored within
%                   the window -- a nonlinear, time-asymmetric variant of
%                   AC1 (cf. SY_RampingWindows.m, which trends this same
%                   statistic across segments rather than summarizing its
%                   spread across windows as done here)
%
% acrossWinStat, controls how the obtained sequence of local estimates is
%                   compared (as a ratio to the full time series):
%                       (i) 'std': standard deviation
%                       (ii) 'ent' (kernel-smoothed) distributional entropy
%                       (iii) 'permen': normalized Permutation Entropy, PermEn(3,1)
%                           of the sequence of local estimates
%
% numSeg, the number of segments to divide the time series up into, thus
%       controlling the window length
%
% incMove, the increment to move the window at each iteration, as 1/fraction of the
%       window length (e.g., incMove = 2, means the window moves half the length of the
%       window at each increment)
%
% NOTE: SY_SlidingWindow(y,'mean','std',X,1) is the same as StatAvX, computed as
%                       SY_StatAv(y,'seg',X);
% cf. "Heart rate control in normal and aborted-SIDS infants", S. M. Pincus et al.
%           Am J. Physiol. Regul. Integr. Comp. Physiol. 264(3) R638 (1993)

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

doPlot = false; % plot outputs

% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(windowStat)
	windowStat = 'mean'; % measure within each window
end
if nargin < 3 || isempty(acrossWinStat)
	acrossWinStat = 'std'; % measure across all windows
end
if nargin < 4 || isempty(numSeg)
	numSeg = 5;
end
if nargin < 5 || isempty(incMove)
	incMove = 2;
end

% ------------------------------------------------------------------------------

winLength = floor(length(y) / numSeg); % length of window
if winLength == 0
	warning('Time-series of length %u is too short for %u windows', length(y), numSeg);
	out = NaN;
	return
end
inc = floor(winLength / incMove); % increment to move at each step
% If increment rounded down to zero, prop it up:
if inc == 0
	inc = 1;
end

numSteps = (floor((length(y) - winLength) / inc) + 1);
qs = zeros(numSteps, 1);

% Convert a step index (stepInd) to a range of indices corresponding to that window:
getWindow = @(stepInd) ((stepInd - 1) * inc + 1:(stepInd - 1) * inc + winLength);

switch windowStat
	case 'mean' % Sliding window mean
		for i = 1:numSteps
			qs(i) = mean(y(getWindow(i)));
		end
	case 'std' % Sliding window std
		for i = 1:numSteps
			qs(i) = std(y(getWindow(i)));
		end
	case 'ent' % Sliding window distributional entropy
		for i = 1:numSteps
			qs(i) = EN_DistributionEntropy(y(getWindow(i)), 'ks', []);
		end
	case 'permen' % Sliding window normalized Permutation Entropy, PermEn(3,1)
		for i = 1:numSteps
			permEn_struct = EN_PermEn(y(getWindow(i)), 3, 1);
			if isstruct(permEn_struct)
				qs(i) = permEn_struct.normPermEn;
			else
				qs(i) = NaN;
			end
		end
	case 'specen' % Sliding window normalized Shannon spectral entropy
		for i = 1:numSteps
			yWin = y(getWindow(i));
			Fy = fft(yWin - mean(yWin));
			Py = abs(Fy(1:floor(end / 2) + 1)).^2; % one-sided power spectrum
			Py = Py(Py > 0); % avoid log(0)
			if length(Py) < 2 || sum(Py) == 0
				qs(i) = NaN;
			else
				Py = Py / sum(Py);
				qs(i) = -sum(Py .* log2(Py)) / log2(length(Py)); % normalized to [0,1]
			end
		end
	case 'mom3' % Third moment
		for i = 1:numSteps
			qs(i) = DN_Moments(y(getWindow(i)), 3, true);
		end
	case 'mom4' % Fourth moment
		for i = 1:numSteps
			qs(i) = DN_Moments(y(getWindow(i)), 4, true);
		end
	case 'mom5' % Fifth moment
		for i = 1:numSteps
			qs(i) = DN_Moments(y(getWindow(i)), 5, true);
		end
	case 'lillie' % Lilliefors test statistic (KS-type distance from a Gaussian CDF)
		% Uses lillietest's raw kstat output directly rather than going via
		% HT_DistributionTest's p-value: MATLAB's lillietest interpolates p from a
		% simulated lookup table and clips to its boundary once the data is clearly
		% non-Gaussian, which saturates the p-value to an identical constant across
		% many real (non-Gaussian) windows. kstat is the underlying continuous
		% statistic the p-value is derived from, so it doesn't have this problem.
		warning('off', 'stats:lillietest:OutOfRangePLow'); warning('off', 'stats:lillietest:OutOfRangePHigh');
		for i = 1:numSteps
			[~, ~, qs(i)] = lillietest(y(getWindow(i)), 0.05, 'norm');
		end
		warning('on', 'stats:lillietest:OutOfRangePLow'); warning('on', 'stats:lillietest:OutOfRangePHigh');
	case 'AC1' % Lag-1 autocorrelation
		for i = 1:numSteps
			qs(i) = CO_AutoCorr(y(getWindow(i)), 1, 'Fourier');
		end
	case 'asymAC1' % Asymmetric, nonlinear AC1 variant: mean(x_t*x_{t+1}^2), z-scored per window
		for i = 1:numSteps
			zw = zscore(y(getWindow(i)));
			qs(i) = mean(zw(1:end - 1) .* zw(2:end).^2);
		end
	otherwise
		error('Unknown statistic ''%s''', windowStat)
end

% Check for all errors (e.g., short time series):
if all(isnan(qs))
	warning('These sliding window settings are not suitable for this time series');
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------------
if doPlot
	figure('color', 'w'); box('on');
	plot(round(winLength / 2):inc:(numSteps - 1) * inc + round(winLength / 2), qs, 'r');
end

% ------------------------------------------------------------------------------
% Compute the output statistic
% ------------------------------------------------------------------------------
switch acrossWinStat
	case 'std'
		% normalized by std of full time series
		out = std(qs) / std(y);
	case 'permen'
		permEn_struct = EN_PermEn(qs, 3, 1); % normalized PermEn of the sliding window measures
		if isstruct(permEn_struct)
			out = permEn_struct.normPermEn;
		else
			out = NaN;
		end
	case 'ent'
		% get a load of statistics from kernel-smoothed distribution (inefficient since only one is used)
		kssimpouts = DN_FitKernelSmooth(qs);
		out = kssimpouts.entropy; % distributional entropy
	otherwise
		error('Unknown statistic: ''%s''.', acrossWinStat)
end

end
