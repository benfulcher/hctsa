function out = SY_DriftingMeanCUSUM(y)
% SY_DriftingMeanCUSUM   Parameter-free CUSUM test for a drifting mean.
%
% SY_StatAv and SY_DriftingMean both test for a drifting mean by splitting the
% time series into segments, requiring a choice of segment length or number of
% segments. This is the analogous test with no free parameter: it forms cumsum(y)
% directly (at full resolution) and computes CUSUM/bridge statistics on it (see
% BF_CumSumBridgeStats), including a comparison between an ordinary least-squares
% and a robust linear fit to flag whether any apparent drift is outlier-driven
% or a genuine trend. The same statistics are applied to lagged products by
% SY_DriftingAutoCorr.
%
% The basic cumsum trend statistics (meanYC, gradient, intercept, meanYC12,
% meanYC22) duplicate what SY_Trend already computes on y's own cumsum (|r| >=
% 0.79 on real EEG data, several essentially exact), so are dropped here.
% stdBridge is also dropped: since the input is always z-scored (mean exactly 0),
% std(bridge) tracks SY_Trend_stdYC almost exactly (r = 1.000 on real-world series).
% Only maxBridge, posMaxBridge, and the robust-vs-OLS comparison statistics,
% which SY_Trend does not provide, are returned.
%
% ---INPUTS:
% y, the input time series (assumed z-scored)
%
% ---OUTPUTS:
% From the cumulative sum yC = cumsum(y), its straight-line 'bridge'
% (yC minus the line joining 0 to its final value), and linear fits to yC; NaN if
% the series is shorter than 20 samples:
% maxBridge, the largest absolute value of the bridge, divided by std(y)*sqrt(N)
% posMaxBridge, the position (from 0 to 1) of that largest value
% gradientDiffSE, the difference between the ordinary least-squares and a robust
%       (bisquare) slope fit to yC, in standard errors of the robust slope
% residStdRatio, the standard deviation of the OLS residuals over that of the
%       robust-fit residuals
% varRatioTrend, the Kendall rank correlation between time and the log of the
%       squared bridge divided by its Brownian-bridge null variance: whether the
%       bridge wanders more (or less) later in the series than expected under
%       stationarity

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

if ~BF_iszscored(y)
	warning('The input time series should be z-scored')
end

out = BF_CumSumBridgeStats(y);
if isstruct(out)
	out = rmfield(out, {'meanYC', 'gradient', 'intercept', 'meanYC12', 'meanYC22', 'stdBridge'});
end

end
