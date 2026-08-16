function out = SY_DriftingAutoCorr(y, tau, whatProduct)
% SY_DriftingAutoCorr  Drift in lag-tau (auto)correlation via a cumulative-sum test.
%
% Forms a lag-tau cross-term product series p_t and tests whether its mean --
% i.e. the corresponding (linear or nonlinear) correlation statistic -- is
% stationary, via CUSUM/bridge statistics on cumsum(p) (see BF_CumSumBridgeStats).
% Under stationarity cumsum(p) grows ~linearly; systematic curvature or a
% localized departure from that line indicates that the correlation structure --
% not just the mean or variance of y itself -- is drifting over the course of
% the time series. This is a CUSUM-style stationarity test (cf. Inclan-Tiao's
% test for a change point in variance) applied to a lag-product series rather
% than to y itself, mirroring what SY_Trend does for y's own cumsum.
%
% ---INPUTS:
% y, the input time series (assumed z-scored)
%
% tau, the lag defining the cross term [default: 1]
%
% whatProduct, which cross term p_t to test for drift:
%              'ac' (default): p_t = y(t).y(t+tau)          -- linear autocorrelation
%              'forward':      p_t = y(t).y(t+tau)^2         -- nonlinear/asymmetric
%                               (does the signed value now predict the squared,
%                               energy-like value later? cf. CO_AutoCorrX2)
%              'backward':     p_t = y(t)^2.y(t+tau)         -- the other direction
%                               (does the squared value now predict the signed
%                               value later?)
%              'asymmetry':    p_t = y(t).y(t+tau).(y(t+tau) - y(t))
%                               = forward - backward, i.e., the leverage/
%                               time-irreversibility signature itself (vanishes
%                               in expectation for time-reversible linear
%                               processes); tests whether *that* asymmetry --
%                               not just its forward or backward half -- is
%                               drifting over the course of the time series.

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

if nargin < 2 || isempty(tau)
	tau = 1;
end
if nargin < 3 || isempty(whatProduct)
	whatProduct = 'ac';
end

if ~BF_iszscored(y)
	warning('The input time series should be z-scored')
end

y = y(:);
N = length(y);

%-------------------------------------------------------------------------------
%% Lag-tau product series
%-------------------------------------------------------------------------------
if tau >= N - 1
	out = NaN;
	return
end
yEarlier = y(1:end - tau); % y(t)
yLater = y(1 + tau:end);   % y(t+tau)
switch whatProduct
	case 'ac'
		p = yEarlier .* yLater;
	case 'forward'
		p = yEarlier .* yLater.^2;
	case 'backward'
		p = yEarlier.^2 .* yLater;
	case 'asymmetry'
		p = yEarlier .* yLater .* (yLater - yEarlier);
	otherwise
		error('Unknown whatProduct ''%s'' (should be ''ac'', ''forward'', ''backward'', or ''asymmetry'')', whatProduct);
end

out = BF_CumSumBridgeStats(p);

end
