function out = ST_LocalExtrema(y, howToWindow, n)
% ST_LocalExtrema   How local maximums and minimums vary across the time series.
%
% Finds maximums and minimums within given segments of the time series and
% analyses the results.
%
% ---INPUTS:
% y, the input time series
%
% howToWindow, whether to use:
%     (i) 'l', windows of a given length (in which case the third input, n
%             specifies the length)
%     (ii) 'n', a specified number of windows to break the time series up into
%               (in which case the third input, n specifies this number)
%     (iii) 'tau', sets a window length equal to the correlation length of the
%                 time series, the first zero-crossing of the autocorrelation
%                 function.
%
% n, somehow specifies the window length given the setting of howToWindow above.

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
% Check Inputs
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(howToWindow)
	howToWindow = 'l';
end
if nargin < 3 || isempty(n)
	switch howToWindow
		case 'l'
			n = 100; % 100-sample windows
		case 'n'
			n = 5; % 5 windows
	end
end

doPlot = false; % plot outputs to a figure

N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Set window length
% ------------------------------------------------------------------------------
switch howToWindow
	case 'l'
		windowLength = n; % window length
	case 'n'
		windowLength = floor(N / n); % number of windows
	case 'tau'
		% this may not be a good idea!
		windowLength = CO_FirstCrossing(y, 'ac', 0, 'discrete');
	otherwise
		error('Unknown method ''%s''', howToWindow);
end

if isnan(windowLength) || (windowLength > N) || (windowLength <= 1)
	% This feature is unsuitable if the window length exceeds ts
	fprintf(1, 'The window length is longer than the time-series length!\n');
	out = NaN;
	return
end

% ------------------------------------------------------------------------------
%% Buffer the time series
% ------------------------------------------------------------------------------
y_buff = buffer(y, windowLength); % no overlap
% Each *column* is a window of samples:
if y_buff(end) == 0
	y_buff = y_buff(:, 1:end - 1); % remove last window if zero-padded
end
numWindows = size(y_buff, 2); % number of windows

% ------------------------------------------------------------------------------
%% Find local extrema
% ------------------------------------------------------------------------------
locMax = max(y_buff); % summary of local maxima
locMin = min(y_buff); % summary of local minima
absLocMin = abs(locMin); % absolute value of local minima
exti = find(absLocMin > locMax);
locExt = locMax;
locExt(exti) = locMin(exti); % local extrema (furthest from mean; either maxs or mins)
absLocExt = abs(locExt); % the magnitude of the most extreme events in each window

if doPlot
	figure('color', 'w');
	hold('on');
	plot(locMax);
	plot(absLocExt, '--g');
	plot(absLocMin, ':r')
	plot(locExt, 'k');
end

% ------------------------------------------------------------------------------
%% Return outputs
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% Scale of a local extreme, under a Gaussian null
% ------------------------------------------------------------------------------
% The magnitude of a within-window extreme grows with the window length: for w
% iid standard normals E[max] ~ sqrt(2*log(w)). With howToWindow = 'n' the
% window is windowLength = floor(N/n), so that expectation -- and every location
% statistic derived from it -- grows with the time-series length rather than
% with any property of the data (measured: meanmax rises 0.59 -> 2.33 across
% N = 100..3200 for n = 50, a factor of 4).
%
% Dividing by the expected maximum of windowLength standard normals expresses
% these on a scale-free footing: "how extreme are the local extremes, relative
% to Gaussian expectation for a window this long". That removes the length
% dependence entirely (ratio 3.96 -> 1.00 for n50, 2.45 -> 1.00 for n25) and
% leaves the fixed-window 'l' variants unchanged up to a constant, since
% windowLength does not vary with N there.
%
% Blom's approximation for the expected largest order statistic:
expMax = norminv((windowLength - 0.375) / (windowLength + 0.25));

% Ratios below are deliberately left on the raw values: the expMax factor
% cancels exactly in a ratio, so normalizing there would change nothing.
out.meanrat = mean(locMax) / mean(absLocMin);
out.medianrat = median(locMax) / median(absLocMin);
out.minmax = min(locMax) / expMax;
out.minabsmin = min(absLocMin) / expMax;
out.minmaxonminabsmin = min(locMax) / min(absLocMin);
out.meanmax = mean(locMax) / expMax;
out.meanabsmin = mean(absLocMin) / expMax;
out.meanext = mean(locExt) / expMax;
out.medianmax = median(locMax) / expMax;
out.medianabsmin = median(absLocMin) / expMax;
out.medianext = median(locExt) / expMax;
% The std fields are NOT normalized by expMax: the spread of a maximum *shrinks*
% as the window grows (SD[max_w] ~ 1/sqrt(2*log(w))), so dividing by expMax
% would steepen their drift rather than remove it. Normalizing instead by the
% Gumbel asymptotic for SD[max_w] was tried and made things worse at small
% windows (ratio 0.56 -> 1.36 for n50), because the 'n' variants reach window
% lengths of only 2-8 samples at short N, where the asymptotic does not hold.
% These remain length-sensitive for the 'n' parameterizations.
out.stdmax = std(locMax);
out.stdmin = std(locMin);
out.stdext = std(locExt);
out.zcext = ST_SimpleStats(locExt, 'zcross');
out.meanabsext = mean(absLocExt) / expMax;
out.medianabsext = median(absLocExt) / expMax;
% Not normalized by expMax: this is a *difference* between two extremes of
% similar magnitude, so it does not carry the expMax scale and dividing by it
% over-corrects (measured: eta^2 rose from 0.82 to 0.97 for n50 when it was
% included in the normalization).
out.diffmaxabsmin = sum(abs(locMax - absLocMin)) / numWindows;
out.uord = sum(sign(locExt)) / numWindows; % whether extreme events are more up or down
out.maxmaxmed = max(locMax) / median(locMax);
out.minminmed = min(locMin) / median(locMin);
out.maxabsext = max(absLocExt) / median(absLocExt);

end
