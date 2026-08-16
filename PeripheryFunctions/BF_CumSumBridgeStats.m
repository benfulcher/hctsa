function out = BF_CumSumBridgeStats(p)
% BF_CumSumBridgeStats  CUSUM/bridge stationarity statistics on a cumulative sum.
%
% Given a series p_t whose mean is being tested for stationarity, forms
% cumsum(p) and computes: linear-fit statistics on the cumsum (cf. SY_Trend),
% a CUSUM 'bridge' relative to the endpoint-to-endpoint line (cf. Inclan-Tiao's
% test for a change point in variance), and a comparison between an ordinary
% least-squares and a robust (bisquare) linear fit to the cumsum -- large
% disagreement between the two indicates the OLS trend is either outlier-driven
% or genuinely curved (accelerating/decelerating drift) rather than a clean
% linear trend.
%
%---INPUT:
% p, the series to test for a stationary mean via its cumulative sum.
%
%---OUTPUT:
% out, a structure of statistics, or NaN if p is too short (< 20 samples).

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

p = p(:);
Np = length(p);

if Np < 20
	out = NaN;
	return
end

t = (1:Np)';
yC = cumsum(p);

out.meanYC = mean(yC);

%-------------------------------------------------------------------------------
%% Ordinary least-squares linear fit to the cumsum
%-------------------------------------------------------------------------------
coeffsOLS = polyfit(t, yC, 1);
out.gradient = coeffsOLS(1); % ~ std(yC) too (r > 0.99 empirically); kept as the interpretable one
out.intercept = coeffsOLS(2);
residOLS = yC - polyval(coeffsOLS, t);

% Mean cumsum in first and second half of the time series (cf. SY_Trend):
out.meanYC12 = mean(yC(1:floor(Np / 2)));
out.meanYC22 = mean(yC(floor(Np / 2) + 1:end));

%-------------------------------------------------------------------------------
%% CUSUM bridge relative to the endpoint-to-endpoint line (Inclan-Tiao-style)
%-------------------------------------------------------------------------------
bridge = yC - (t / Np) * yC(end);
scaleFactor = std(p) * sqrt(Np);
if scaleFactor > 0
	out.maxBridge = max(abs(bridge)) / scaleFactor;
else
	out.maxBridge = NaN;
end
[~, idxMax] = max(abs(bridge));
out.posMaxBridge = idxMax / Np; % where the largest deviation from stationarity occurs
out.stdBridge = std(bridge);

%-------------------------------------------------------------------------------
%% Robust vs. OLS regression: is the OLS trend outlier-driven or genuine drift?
%-------------------------------------------------------------------------------
warning('off', 'stats:statrobustfit:IterationLimit')
[robCoeffs, robStats] = robustfit(t, yC);
warning('on', 'stats:statrobustfit:IterationLimit')
robustGradient = robCoeffs(2); % not output directly: r > 0.98 with out.gradient
if robStats.se(2) > 0
	out.gradientDiffSE = (out.gradient - robustGradient) / robStats.se(2);
else
	out.gradientDiffSE = NaN;
end
if std(robStats.resid) > 0
	out.residStdRatio = std(residOLS) / std(robStats.resid);
else
	out.residStdRatio = NaN;
end

end
