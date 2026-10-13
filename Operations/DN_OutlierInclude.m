function out = DN_OutlierInclude(y, thresholdHow, inc, fixedThresh)
% DN_OutlierInclude   How the timing and spacing of extreme values change as the threshold rises.
%
% Raises a threshold th from 0 to the maximum value of the series, in increments
% of inc, and at each threshold takes the "events": the points at or beyond it
% (for 'abs', values with abs(y) >= th; for 'pos', y >= th; for 'neg',
% y <= -th). The threshold is applied to y itself, so the series should be
% z-scored (a warning is issued otherwise). At each threshold it records:
%   (1) the mean gap (in samples) between successive events, and its standard
%       error (std of the gaps / sqrt of their number),
%   (2) the percentage of points that are events (the number of gaps over the
%       number of candidate points, times 100),
%   (3) the median and mean time of the events, rescaled so that the start of the
%       series is -1, the middle is 0 and the end is 1, and std(times)/sqrt(their
%       number) (in samples).
% The sweep stops when events are 2% or fewer of the points. The outputs measure
% how these curves change with th, using exponential [f(x) = a*exp(b*x) + c] and
% linear [f(x) = a*x + b] least-squares fits, and simple statistics across
% thresholds. The exponential fit is the global optimum over the rate b, with a and
% c found by linear least squares (see BF_ExpFit), so it does not depend on a
% starting point or optimizer settings; its R^2 lies between 0 and 1. If a fit is
% not meaningful (too few thresholds, or a constant curve), the corresponding
% outputs are NaN.
%
% If fixedThresh is given, the sweep and fits are skipped, and the statistics in
% (1)-(3) are returned for that one threshold.
%
% ---INPUTS:
% y, the input time series (ideally z-scored)
%
% thresholdHow, how to determine outliers:
%       (i) 'abs': values furthest from zero in either direction (default),
%       (ii) 'pos': the greatest positive values, or
%       (iii) 'neg': the greatest negative values.
%
% inc, the increment to move through (in units of the standard deviation if the
%       time series is z-scored; default 0.01). Unused when fixedThresh is given.
%
% fixedThresh, [optional] a single threshold (in the units of y, e.g., 2 for two
%       standard deviations of a z-scored series). If given, the sweep is skipped.
%
% ---OUTPUTS:
% From the sweep (fixedThresh not given):
% mfexpb, mfexpr2, mfexprmse: the rate b, R^2 and root-mean-square error of the
%       exponential fit to the mean gap vs. th (the amplitude a and offset c are not
%       output: they become large and opposite in sign when the curve is nearly
%       straight)
% nfexpb, nfexpr2, nfexprmse: the same for an exponential fit to the percentage of
%       points that are events vs. th
% nfla, nflb, nflr2, nflrmse: slope a, intercept b, R^2 and root-mean-square
%       error of a linear fit to the percentage of points that are events vs. th
% mdtm, mdtmd, mdtstd: mean, median and standard deviation of the mean gap
%       across thresholds
% mdrm, mdrmd, mdrstd: mean, median and standard deviation, across thresholds,
%       of the median time of the events (-1 to 1)
% mrm, mrmd, mrstd: mean, median and standard deviation, across thresholds, of
%       the mean time of the events (-1 to 1)
% xcmerr1, xcmerrn1: cross-correlation (xcorr with 'coeff' normalization)
%       between the mean gap and its standard error across thresholds, at lags
%       +1 and -1
% stdrfexpb, stdrfexpr2, stdrfexprmse: the rate and fit quality of an exponential
%       fit to std(times)/sqrt(their number) vs. th
% stdrfla, stdrflb, stdrflr2, stdrflrmse: the same for a linear fit
% From a single threshold (fixedThresh given; all NaN except propIncluded if events
% are 2% or fewer of the points):
% meanDt, seDt: the mean gap between events, and its standard error
% propIncluded: the percentage of points that are events
% medianRelTime, meanRelTime: the median and mean time of the events (-1 to 1)
% stdRelTime: std(times)/sqrt(their number), in samples
%
% A constant time series returns NaN.

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
%% Preliminaries
% ------------------------------------------------------------------------------

doPlot = false; % Plot some outputs

% ------------------------------------------------------------------------------
%% Check inputs and set defaults
% ------------------------------------------------------------------------------
% If the time series is a constant causes issues
if all(y == y(1))
	% This method is not suitable for such time series: return a NaN
	fprintf(1, 'The time series is a constant!\n');
	out = NaN;
	return
end

% Check z-scored time series
if ~BF_iszscored(y)
	warning('The input time series should be z-scored')
end
N = length(y); % length of the time series

if nargin < 2 || isempty(thresholdHow)
	thresholdHow = 'abs'; % Analyze absolute value deviations in the time series by default
end

if nargin < 4
	fixedThresh = [];
end

% ------------------------------------------------------------------------------
%% Fixed-threshold fast path
% ------------------------------------------------------------------------------
% Bypasses the sweep-and-curve-fit machinery below entirely: no Curve Fitting
% toolbox needed, no NLS convergence to worry about. Evaluates the same
% exceedance/interval statistics the sweep loop computes at each threshold,
% but at just this one threshold, and returns them directly rather than
% fitting a trend to how they change across many thresholds.
if ~isempty(fixedThresh)
	switch thresholdHow
		case 'abs'
			r = find(abs(y) >= fixedThresh);
			tot = N;
		case 'pos'
			r = find(y >= fixedThresh);
			tot = sum(y >= 0);
		case 'neg'
			r = find(y <= -fixedThresh);
			tot = sum(y <= 0);
		otherwise
			error('Error thresholding with ''%s''. Must select either ''abs'', ''pos'', or ''neg''.', thresholdHow)
	end

	Dt_exc = diff(r);
	meanDt = mean(Dt_exc);
	propIncluded = length(Dt_exc) / tot * 100;

	% Same "too few events to say anything meaningful" bar as the sweep loop
	% (trimthr = 2%) below, so results from the two modes stay comparable:
	if isnan(meanDt) || propIncluded <= 2
		out.meanDt = NaN;
		out.seDt = NaN;
		out.propIncluded = propIncluded;
		out.medianRelTime = NaN;
		out.meanRelTime = NaN;
		out.stdRelTime = NaN;
		return
	end

	out.meanDt = meanDt;
	out.seDt = std(Dt_exc) / sqrt(length(Dt_exc));
	out.propIncluded = propIncluded;
	out.medianRelTime = median(r) / (N / 2) - 1; % median event timing, relative to middle (N/2)
	out.meanRelTime = mean(r) / (N / 2) - 1; % mean event timing, relative to middle (N/2)
	out.stdRelTime = std(r) / sqrt(length(r));
	return
end

if nargin < 3 || isempty(inc)
	inc = 0.01; % increment through z-scored time-series values
end

% ------------------------------------------------------------------------------
%% Initialize thresholds
% ------------------------------------------------------------------------------
% Could be better to just use a fixed number of increments here, from 0 to the max.
% (rather than forcing a fixed inc)
switch thresholdHow
	case 'abs' % analyze absolute value deviations
		thr = (0:inc:max(abs(y)));
		tot = N;
	case 'pos' % analyze only positive deviations
		thr = (0:inc:max(y));
		tot = sum(y >= 0);
	case 'neg' % analyze only negative deviations
		thr = (0:inc:max(-y));
		tot = sum(y <= 0);
	otherwise
		error('Error thresholding with ''%s''. Must select either ''abs'', ''pos'', or ''neg''.', thresholdHow)
end

if isempty(thr)
	error('Error setting increments through the time-series values...')
end

% -------------------------------------------------------------------------------
% Calculate statistics of over-threshold events, looping over thresholds
% -------------------------------------------------------------------------------
% Break out of the loop as soon as too few points remain to be useful,
% instead of computing statistics across the full threshold range and
% trimming afterward: raising the threshold can only shrink the qualifying
% set of points (never grow it), so once the "at least 2% of data included"
% criterion fails (or too few points remain to form any interval), it will
% keep failing for every higher threshold too -- continuing past that point
% is pure wasted work.
%
% Also searches only within the previous threshold's qualifying indices
% (rPrev) instead of rescanning the full series y at every threshold: since
% the qualifying set only shrinks as the threshold rises, any index that
% didn't qualify at a lower threshold can't qualify at a higher one either.
% This also preserves ascending time-index order for free (logical indexing
% preserves relative order), which matters because Dt_exc = diff(r) and the
% median(r)/mean(r) stats below depend on r being in temporal order, not
% sorted by value.
% -------------------------------------------------------------------------------
trimthr = 2; % percent

msDt = zeros(length(thr), 6); % mean, std, proportion_of_time_series_included,
% median of index relative to middle, mean,
% error
rPrev = (1:N)'; % candidate index set, narrows each iteration
numValid = 0;
for i = 1:length(thr)
	th = thr(i); % the threshold

	% Construct a series consisting of inter-event intervals for parts
	% of the time series exceeding the threshold, th (in a given direction).
	% Searches only the previous iteration's qualifying indices, rPrev:
	switch thresholdHow
		case 'abs' % look at absolute value deviations
			r = rPrev(abs(y(rPrev)) >= th);
		case 'pos' % look at only positive deviations
			r = rPrev(y(rPrev) >= th);
		case 'neg' % look at only negative deviations
			r = rPrev(y(rPrev) <= -th);
	end
	rPrev = r; % narrower candidate set for the next (higher) threshold

	% The Dt (interval) time series of values exceeding threshold
	Dt_exc = diff(r);

	% ~~~~~~~~~~~~
	% Statistics on the interval time series:
	% ~~~~~~~~~~~~
	meanDt = mean(Dt_exc); % the mean value of inter-event intervals
	propIncluded = length(Dt_exc) / tot * 100; % this is just really measuring the distribution
	% : the proportion of possible values
	% that are actually used in
	% calculation

	% Stop once too few events remain to say anything meaningful (Dt_exc
	% empty/singleton -> NaN mean), or once statistical power is lacking
	% (less than 2% of data included):
	if isnan(meanDt) || propIncluded <= trimthr
		break
	end

	numValid = numValid + 1;
	msDt(numValid, 1) = meanDt;
	msDt(numValid, 2) = std(Dt_exc) / sqrt(length(Dt_exc)); % error on the mean
	msDt(numValid, 3) = propIncluded;
	% ~~~~~~~~~~~~
	% Statistics on the indices of over-threshold events:
	% ~~~~~~~~~~~~
	% The [x/(N/2)-1] rescales such that the middle index, N/2 => 0, and N maps to 1, 0 maps to -1.
	msDt(numValid, 4) = median(r) / (N / 2) - 1; % the median timing of events (relative to middle, N/2)
	msDt(numValid, 5) = mean(r) / (N / 2) - 1; % mean timing of events (relative to middle, N/2)
	msDt(numValid, 6) = std(r) / sqrt(length(r)); % variance of event timing
end
msDt = msDt(1:numValid, :);
thr = thr(1:numValid);

% ------------------------------------------------------------------------------
%% Plot output
% ------------------------------------------------------------------------------
if doPlot
	figure('color', 'w');
	hold('on')
	plot(thr, msDt(:, 1), '.-k');
	plot(thr, msDt(:, 2), '.-b');
	plot(thr, msDt(:, 3), '.-g');
	plot(thr, msDt(:, 4) * 100, '.-m');
	plot(thr, msDt(:, 5) * 100, '.-r');
	plot(thr, msDt(:, 6), '.-c'); hold off
end

% -------------------------------------------------------------------------------
% Quantify outputs:
% -------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Fit an exponential to the mean inter-event interval as a function of the threshold
% ------------------------------------------------------------------------------
% (Global least-squares fits by variable projection: see BF_ExpFit. The amplitude a
% and offset c are not output: they are poorly determined when the curve is close
% to a straight line, as a and c then become large and opposite in sign)
fExp = BF_ExpFit(thr', msDt(:, 1), true);
out.mfexpb = fExp.b;
out.mfexpr2 = fExp.r2;
out.mfexprmse = fExp.rmse;

% ------------------------------------------------------------------------------
%% Fit an exponential to N: the valid proportion left in calculation
% ------------------------------------------------------------------------------
fExp = BF_ExpFit(thr', msDt(:, 3), true);
out.nfexpb = fExp.b;
out.nfexpr2 = fExp.r2;
out.nfexprmse = fExp.rmse;

% ------------------------------------------------------------------------------
%% Fit an linear trend to N: the valid proportion left in calculation
% ------------------------------------------------------------------------------
fLin = SUB_LinearFit(thr', msDt(:, 3));
out.nfla = fLin.a;
out.nflb = fLin.b;
out.nflr2 = fLin.r2;
out.nflrmse = fLin.rmse;

% ------------------------------------------------------------------------------
%% Stationarity metrics
% ------------------------------------------------------------------------------
% Mean, median and std of the mean inter-event interval:
out.mdtm = mean(msDt(:, 1));
out.mdtmd = median(msDt(:, 1));
out.mdtstd = std(msDt(:, 1));

% Mean, median and std of the median and mean of indices over-threshold events occur
out.mdrm = mean(msDt(:, 4));
out.mdrmd = median(msDt(:, 4));
out.mdrstd = std(msDt(:, 4));

out.mrm = mean(msDt(:, 5));
out.mrmd = median(msDt(:, 5));
out.mrstd = std(msDt(:, 5));

% ------------------------------------------------------------------------------
%% Cross correlation between mean and error
% ------------------------------------------------------------------------------
xc = xcorr(msDt(:, 1), msDt(:, 2), 1, 'coeff');
out.xcmerr1 = xc(end); % this is the cross-correlation at lag 1
out.xcmerrn1 = xc(1); % this is the cross-correlation at lag -1

% ------------------------------------------------------------------------------
%% Fit exponential to std in range
% ------------------------------------------------------------------------------
fExp = BF_ExpFit(thr', msDt(:, 6), true);
out.stdrfexpb = fExp.b;
out.stdrfexpr2 = fExp.r2;
out.stdrfexprmse = fExp.rmse;

% ------------------------------------------------------------------------------
%% Fit linear to errors in range
% ------------------------------------------------------------------------------
fLin = SUB_LinearFit(thr', msDt(:, 6));
out.stdrfla = fLin.a;
out.stdrflb = fLin.b;
out.stdrflr2 = fLin.r2;
out.stdrflrmse = fLin.rmse;

if doPlot
	figure('color', 'w')
	errorbar(thr, msDt(:, 1), msDt(:, 2), 'k');
end

end

% ------------------------------------------------------------------------------
function f = SUB_LinearFit(x, y)
	% Ordinary least-squares line y = a*x + b, with R^2 and root-mean-square error (n - 2 degrees of freedom)
	n = length(y);
	if n < 3 || ~(sum((y - mean(y)).^2) > 0)
		% Too few points, or a constant curve: no meaningful fit
		f = struct('a', NaN, 'b', NaN, 'r2', NaN, 'rmse', NaN);
		return
	end
	p = polyfit(x, y, 1);
	SSE = sum((y - polyval(p, x)).^2);
	SST = sum((y - mean(y)).^2);
	f.a = p(1);
	f.b = p(2);
	f.r2 = 1 - SSE/SST;
	f.rmse = sqrt(SSE/(n - 2));
end
