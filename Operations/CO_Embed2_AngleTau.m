function out = CO_Embed2_AngleTau(y, maxTau)
% CO_Embed2_AngleTau   Angle autocorrelation in a 2-dimensional delay embedding, over delays.
%
% Investigates how the autocorrelation of the angles between successive points in
% the two-dimensional delay embedding (y(t), y(t+tau)) changes as tau varies from
% 1, 2, ..., maxTau. For each tau the angle of each step is measured from the
% horizontal, as in CO_Embed2, and its autocorrelation is taken at lags 1, 2 and 3.
% The outputs summarize the three resulting curves as a function of tau.
%
% ---INPUTS:
% y, a column vector time series
% maxTau, the maximum time delay to consider
%
% ---OUTPUTS:
% ac1_thetaac1, ac1_thetaac2, ac1_thetaac3, the lag-1 autocorrelation (across tau)
%       of the lag-1, lag-2 and lag-3 angle autocorrelations,
% mean_thetaac1, mean_thetaac2, mean_thetaac3, their means across tau,
% max_thetaac1, max_thetaac2, max_thetaac3, their maxima,
% min_thetaac1, min_thetaac2, min_thetaac3, their minima,
% meanrat_thetaac12, the ratio mean_thetaac1/mean_thetaac2,
% diff_thetaac12, the sum across tau of the absolute difference between the lag-2
%       and lag-1 angle autocorrelations.
% The output is a single NaN if the series is too short to embed.

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

doPlot = false;
tauRange = (1:1:maxTau);
numTau = length(tauRange);

% Ensure y is a column vector
if size(y, 2) > size(y, 1);
	y = y';
end

stats_store = zeros(3, numTau);

for i = 1:numTau
	tau = tauRange(i);

	m = [y(1:end - tau), y(1 + tau:end)];

	theta = diff(m(:, 2)) ./ diff(m(:, 1));
	theta = atan(theta); % measured as deviation from the horizontal

	if isempty(theta)
		% Data-dependent (too short to embed at this tau), so NaN, not error():
		warning('Time series (N=%u) too short for embedding at tau = %u', length(y), tau);
		out = NaN; return
	end

	stats_store(1, i) = CO_AutoCorr(theta, 1, 'Fourier');
	stats_store(2, i) = CO_AutoCorr(theta, 2, 'Fourier');
	stats_store(3, i) = CO_AutoCorr(theta, 3, 'Fourier');
end

if doPlot
	figure('color', 'w'); box('on');
	plot(stats_store');
end

% ------------------------------------------------------------------------------
% Compute lots of outputs statistics:
% ------------------------------------------------------------------------------
out.ac1_thetaac1 = CO_AutoCorr(stats_store(1, :), 1, 'Fourier');
out.ac1_thetaac2 = CO_AutoCorr(stats_store(2, :), 1, 'Fourier');
out.ac1_thetaac3 = CO_AutoCorr(stats_store(3, :), 1, 'Fourier');
out.mean_thetaac1 = mean(stats_store(1, :));
out.max_thetaac1 = max(stats_store(1, :));
out.min_thetaac1 = min(stats_store(1, :));
out.mean_thetaac2 = mean(stats_store(2, :));
out.max_thetaac2 = max(stats_store(2, :));
out.min_thetaac2 = min(stats_store(2, :));
out.mean_thetaac3 = mean(stats_store(3, :));
out.max_thetaac3 = max(stats_store(3, :));
out.min_thetaac3 = min(stats_store(3, :));
out.meanrat_thetaac12 = out.mean_thetaac1 / out.mean_thetaac2;
out.diff_thetaac12 = sum(abs(stats_store(2, :) - stats_store(1, :)));

end
