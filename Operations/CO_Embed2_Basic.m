function out = CO_Embed2_Basic(y, tau)
% CO_Embed2_Basic   Point-density statistics in a 2-d embedding space.
%
% Plots y(t+tau) against y(t) and computes the fraction of points that fall inside
% simple geometric shapes: thick diagonals, thick parabolas, rings and circles. The
% input is assumed to be z-scored.
%
% ---INPUTS:
% y, the input time series
% tau, the time lag (default 1): a number of samples, or a string that sets it from
%      the series: 'tau' (the first zero-crossing of the autocorrelation function,
%      capped at N/10), 'ac1e' (the floor of its first 1/e crossing), or 'mi' (the
%      smaller of the first minimum of the Kraskov automutual information and the
%      'ac1e' delay); see BF_GetTau
%
% ---OUTPUTS:
% updiag01, updiag05, the fraction of points within 0.1 or 0.5 (vertically) of the
%       diagonal y(t+tau) = y(t),
% downdiag01, downdiag05, the same for the anti-diagonal y(t+tau) = -y(t),
% ratdiag01, ratdiag05, the ratios updiag/downdiag,
% parabup01, parabup05, parabdown01, parabdown05, the fraction within 0.1 or 0.5 of
%       the parabolas y(t+tau) = y(t)^2 and y(t+tau) = -y(t)^2,
% parabup01_1, parabup05_1, parabdown01_1, parabdown05_1, the same for the
%       parabolas shifted up by 1 (y(t)^2 + 1 and -(y(t)^2 - 1)),
% parabup01_n1, parabup05_n1, parabdown01_n1, parabdown05_n1, the same for the
%       parabolas shifted down by 1 (y(t)^2 - 1 and -(y(t)^2 + 1)),
% ring1_01, ring1_02, ring1_05, the fraction with |y(t)^2 + y(t+tau)^2 - 1| below
%       0.1, 0.2 or 0.5,
% incircle_01, incircle_02, incircle_05, incircle_1, incircle_2, incircle_3, the
%       fraction with y(t)^2 + y(t+tau)^2 below 0.1, 0.2, 0.5, 1, 2 or 3,
% medianincircle, stdincircle, the median and standard deviation of the six incircle
%       fractions.

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

if nargin < 2
	tau = 1;
end

doPlot = false; % plot outputs to a figure

if ischar(tau) && ismember(tau, {'ac1e', 'mi'})
	% Adaptive delay: see BF_GetTau
	tau = BF_GetTau(y, tau);
	if isnan(tau)
		out = NaN; return
	end
end
if strcmp(tau, 'tau')
	% Make tau the first zero crossing of the autocorrelation function
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
	if isnan(tau)
		out = NaN; return
	end
	% Cannot set the time delay greater than 10% the length of the time series
	if tau > length(y) / 10
		tau = floor(length(y) / 10);
	end
end
if isnan(tau)
	out = NaN; return
end

xt = y(1:end - tau); % part of the time series
xtp = y(1 + tau:end); % time-lagged time series
N = length(y) - tau; % Length of each time series subsegment

% Points in a thick bottom-left -- top-right diagonal
out.updiag01 = sum(abs(xtp - xt) < 0.1) / N;
out.updiag05 = sum(abs(xtp - xt) < 0.5) / N;

% Points in a thick bottom-right -- top-left diagonal
out.downdiag01 = sum(abs(xtp + xt) < 0.1) / N;
out.downdiag05 = sum(abs(xtp + xt) < 0.5) / N;

% Ratio of these
out.ratdiag01 = out.updiag01 / out.downdiag01;
out.ratdiag05 = out.updiag05 / out.downdiag05;

% In a thick parabola concave up
out.parabup01 = sum(abs(xtp - xt.^2) < 0.1) / N;
out.parabup05 = sum(abs(xtp - xt.^2) < 0.5) / N;

% In a thick parabola concave down
out.parabdown01 = sum(abs(xtp + xt.^2) < 0.1) / N;
out.parabdown05 = sum(abs(xtp + xt.^2) < 0.5) / N;

% In a thick parabola concave up, shifted up 1
out.parabup01_1 = sum(abs(xtp - (xt.^2 + 1)) < 0.1) / N;
out.parabup05_1 = sum(abs(xtp - (xt.^2 + 1)) < 0.5) / N;

% In a thick parabola concave down, shifted up 1
out.parabdown01_1 = sum(abs(xtp + (xt.^2 - 1)) < 0.1) / N;
out.parabdown05_1 = sum(abs(xtp + (xt.^2 - 1)) < 0.5) / N;

% In a thick parabola concave up, shifted down 1
out.parabup01_n1 = sum(abs(xtp - (xt.^2 - 1)) < 0.1) / N;
out.parabup05_n1 = sum(abs(xtp - (xt.^2 - 1)) < 0.5) / N;

% In a thick parabola concave down, shifted down 1
out.parabdown01_n1 = sum(abs(xtp + (xt.^2 + 1)) < 0.1) / N;
out.parabdown05_n1 = sum(abs(xtp + (xt.^2 + 1)) < 0.5) / N;

% RINGS (points within a radius range)
out.ring1_01 = sum(abs(xtp.^2 + xt.^2 - 1) < 0.1) / N;
out.ring1_02 = sum(abs(xtp.^2 + xt.^2 - 1) < 0.2) / N;
out.ring1_05 = sum(abs(xtp.^2 + xt.^2 - 1) < 0.5) / N;

% CIRCLES (points inside a given circular boundary)
out.incircle_01 = sum(xtp.^2 + xt.^2 < 0.1) / N;
out.incircle_02 = sum(xtp.^2 + xt.^2 < 0.2) / N;
out.incircle_05 = sum(xtp.^2 + xt.^2 < 0.5) / N;
out.incircle_1 = sum(xtp.^2 + xt.^2 < 1) / N;
out.incircle_2 = sum(xtp.^2 + xt.^2 < 2) / N;
out.incircle_3 = sum(xtp.^2 + xt.^2 < 3) / N;
out.medianincircle = median([out.incircle_01, out.incircle_02, out.incircle_05 ...
							 out.incircle_1, out.incircle_2, out.incircle_3]);
out.stdincircle = std([out.incircle_01, out.incircle_02, out.incircle_05 ...
					   out.incircle_1, out.incircle_2, out.incircle_3]);

if doPlot
	figure('color', 'w'); box('on');
	plot(xt, xtp, '.k');
	hold on
	r = (xtp.^2 + xt.^2 < 0.2);
	plot(xt(r), xtp(r), '.g')
end

end
