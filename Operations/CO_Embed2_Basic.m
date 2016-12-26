function out = CO_Embed2_Basic(y,tau)
% CO_Embed2_Basic
%
% Computes a set of point density measures in a plot of y_i against y_{i-tau}.
%
% INPUTS:
% y, the input time series
%
% tau, the time lag (can be set to 'tau' to set the time lag the first zero
%                       crossing of the autocorrelation function)
%
% Outputs include the number of points near the diagonal, and similarly, the
% number of points that are close to certain geometric shapes in the y_{i-tau},
% y_{tau} plot, including parabolas, rings, and circles.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

doPlot = 0; % plot outputs to a figure

if strcmp(tau,'tau')
	% Make tau the first zero crossing of the autocorrelation function
    tau = CO_FirstZero(y,'ac');
end

xt = y(1:end-tau); % part of the time series
xtp = y(1+tau:end); % time-lagged time series
N = length(y) - tau; % Length of each time series subsegment

% Points in a thick bottom-left -- top-right diagonal
out.updiag01 = sum(abs(xtp-xt) < 0.1)/N;
out.updiag05 = sum(abs(xtp-xt) < 0.5)/N;

% Points in a thick bottom-right -- top-left diagonal
out.downdiag01 = sum(abs(xtp+xt) < 0.1)/N;
out.downdiag05 = sum(abs(xtp+xt) < 0.5)/N;

% Ratio of these
out.ratdiag01 = out.updiag01/out.downdiag01;
out.ratdiag05 = out.updiag05/out.downdiag05;

% In a thick parabola concave up
out.parabup01 = sum(abs(xtp-xt.^2) < 0.1)/N;
out.parabup05 = sum(abs(xtp-xt.^2) < 0.5)/N;

% In a thick parabola concave down
out.parabdown01 = sum(abs(xtp+xt.^2) < 0.1)/N;
out.parabdown05 = sum(abs(xtp+xt.^2) < 0.5)/N;

% In a thick parabola concave up, shifted up 1
out.parabup01_1 = sum(abs(xtp-(xt.^2+1)) < 0.1)/N;
out.parabup05_1 = sum(abs(xtp-(xt.^2+1)) < 0.5)/N;

% In a thick parabola concave down, shifted up 1
out.parabdown01_1 = sum(abs(xtp+(xt.^2-1)) < 0.1)/N;
out.parabdown05_1 = sum(abs(xtp+(xt.^2-1)) < 0.5)/N;

% In a thick parabola concave up, shifted down 1
out.parabup01_n1 = sum(abs(xtp-(xt.^2-1)) < 0.1)/N;
out.parabup05_n1 = sum(abs(xtp-(xt.^2-1)) < 0.5)/N;

% In a thick parabola concave down, shifted down 1
out.parabdown01_n1 = sum(abs(xtp+(xt.^2+1)) < 0.1)/N;
out.parabdown05_n1 = sum(abs(xtp+(xt.^2+1)) < 0.5)/N;

% RINGS (points within a radius range)
out.ring1_01 = sum(abs(xtp.^2+xt.^2-1) < 0.1)/N;
out.ring1_02 = sum(abs(xtp.^2+xt.^2-1) < 0.2)/N;
out.ring1_05 = sum(abs(xtp.^2+xt.^2-1) < 0.5)/N;

% CIRCLES (points inside a given circular boundary)
out.incircle_01 = sum(xtp.^2+xt.^2 < 0.1)/N;
out.incircle_02 = sum(xtp.^2+xt.^2 < 0.2)/N;
out.incircle_05 = sum(xtp.^2+xt.^2 < 0.5)/N;
out.incircle_1 = sum(xtp.^2+xt.^2 < 1)/N;
out.incircle_2 = sum(xtp.^2+xt.^2 < 2)/N;
out.incircle_3 = sum(xtp.^2+xt.^2 < 3)/N;
out.medianincircle = median([out.incircle_01, out.incircle_02, out.incircle_05 ...
                            out.incircle_1, out.incircle_2, out.incircle_3]);
out.stdincircle = std([out.incircle_01, out.incircle_02, out.incircle_05 ...
                        out.incircle_1, out.incircle_2, out.incircle_3]);

if doPlot
    figure('color','w'); box('on');
    plot(xt,xtp,'.k');
    hold on
    r = (xtp.^2+xt.^2 < 0.2);
    plot(xt(r),xtp(r),'.g')
end

end
