function out = WL_fBM(y)
% WL_fBM   Hurst exponent of the integrated time series treated as a fractional Brownian motion.
%
% Uses the wfbmesti function from Matlab's Wavelet Toolbox, which assumes the input
% is a path of fractional Brownian motion (fBm). The series is treated as the
% increments of that path (fractional Gaussian noise, fGn), so its cumulative sum
% is what is passed to wfbmesti, and H is the Hurst exponent of the fGn: 0.5 for
% white noise, above 0.5 for persistent series, below 0.5 for anti-persistent ones.
% A series that is itself a path (e.g., a random walk), or is strongly
% autocorrelated, has H above 1 (the estimator then gives about 1.2-1.3), outside
% the range of fractional Brownian motion.
%
% ---INPUTS:
% y, the time series to analyze
%
% ---OUTPUTS:
% H_deriv2, the Hurst exponent from wfbmesti's second-order discrete derivative
%        estimator
% H_deriv2Wavelet, the Hurst exponent from the wavelet-based version of the same
%        estimator (using a fixed sym5 filter)
%
% ---NOTES:
% wfbmesti's third estimator (a wavelet-variance-vs-level regression, fixed to a Haar
% decomposition regardless of the wavelet used elsewhere in this codebase) is not
% returned here: WL_modwtvar's decaySlope estimates the same quantity via the
% MODWT's unbiased, boundary-corrected variance decomposition, a more principled
% route to the same scaling exponent.

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
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox');

% Parameter estimation of fractional Brownian motion
% wfbmesti expects an fBm path, not its increments. Simulated fGn with known H
% (N = 1000, 30 series each at H = 0.2, 0.35, 0.5, 0.65, 0.8, 0.95) gives
% H_deriv2 of -0.24, -0.11, 0.01, 0.12, 0.23, 0.34 when passed directly, but
% 0.19, 0.35, 0.51, 0.66, 0.80, 0.94 when its cumulative sum is passed.
hest = wfbmesti(cumsum(y));
out.H_deriv2 = hest(1); % second-order discrete-derivative estimate
out.H_deriv2Wavelet = hest(2); % second-order discrete derivative, wavelet (sym5) version

end
