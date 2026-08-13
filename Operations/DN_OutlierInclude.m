function out = DN_OutlierInclude(y, thresholdHow, inc, fixedThresh)
% DN_OutlierInclude     How statistics depend on distributional outliers.
%
% Measures a range of different statistics about the time series as more and
% more outliers are included in the calculation according to a specified rule,
% of outliers being furthest from the mean, greatest positive, or negative
% deviations.
%
% The threshold for including time-series data points in the analysis increases
% from zero to the maximum deviation, in increments of 0.01*sigma (by default),
% where sigma is the standard deviation of the time series.
%
% At each threshold, the mean, standard error, proportion of time series points
% included, median, and standard deviation are calculated, and outputs from the
% algorithm measure how these statistical quantities change as more extreme
% points are included in the calculation.
%
% ---INPUTS:
% y, the input time series (ideally z-scored)
%
% thresholdHow, the method of how to determine outliers:
%     (i) 'abs': outliers are furthest from the mean,
%     (ii) 'pos': outliers are the greatest positive deviations from the mean, or
%     (iii) 'neg': outliers are the greatest negative deviations from the mean.
%
% inc, the increment to move through (fraction of std if input time series is
%       z-scored). Unused when fixedThresh is given.
%
% fixedThresh [opt], if given, skips the threshold sweep entirely and
%       instead evaluates the same event/interval statistics at this single
%       threshold (a scalar, in the same units as inc, e.g. 2 for two
%       standard deviations of a z-scored series). Returns a small,
%       curve-fit-free struct (meanDt, seDt, propIncluded, and the
%       median/mean/std of event timing) rather than the swept exponential/
%       linear trend fits below, which measure something different: how
%       those quantities *change* as the threshold rises, not their value
%       at one threshold.
%
% Most of the outputs measure either exponential [f(x) = Aexp(Bx) + C] or
% linear [f(x) = Ax + B] fits to the sequence of statistics obtained in
% this way.
%
% [future: could compare differences in outputs obtained with 'p', 'n', and
%               'abs' -- could give an idea as to asymmetries/nonstationarities??]

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

% Check a Curve Fitting toolbox license is available (only needed for the
% sweep-and-fit path below; the fixed-threshold path above has already
% returned by this point):
BF_CheckToolbox('curve_fitting_toolbox');

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
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [0.1 2.5 1]);
f = fittype('a*exp(b*x)+c', 'options', s);
emsg = '';
try
	[c, gof] = fit(thr', msDt(:, 1), f);
catch emsg
	fprintf(1, 'DN_OutlierInclude: error fitting exponential growth to means: %s\n', emsg);
end

if isempty(emsg)
	out.mfexpa = c.a;
	out.mfexpb = c.b;
	out.mfexpc = c.c;
	out.mfexpr2 = gof.rsquare;
	out.mfexprmse = gof.rmse;
else
	out.mfexpa = NaN;
	out.mfexpb = NaN;
	out.mfexpc = NaN;
	out.mfexpr2 = NaN;
	out.mfexprmse = NaN;
end

% ------------------------------------------------------------------------------
%% Fit an exponential to N: the valid proportion left in calculation
% ------------------------------------------------------------------------------
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [120, -1, -16]);
f = fittype('a*exp(b*x)+c', 'options', s);
emsg = '';
try
	[c, gof] = fit(thr', msDt(:, 3), f);
catch emsg
	fprintf(1, 'DN_OutlierInclude: error fitting exponential decay to valid proportion: %s\n', emsg);
end

if isempty(emsg)
	out.nfexpa = c.a;
	out.nfexpb = c.b;
	out.nfexpc = c.c; % (is linearly anticorrelated with c.a)
	out.nfexpr2 = gof.rsquare;
	out.nfexprmse = gof.rmse;
else
	out.nfexpa = NaN;
	out.nfexpb = NaN;
	out.nfexpc = NaN;
	out.nfexpr2 = NaN;
	out.nfexprmse = NaN;
end

% ------------------------------------------------------------------------------
%% Fit an linear trend to N: the valid proportion left in calculation
% ------------------------------------------------------------------------------
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [-40, 100]);
f = fittype('a*x+b', 'options', s);
emsg = '';
try
	[c, gof] = fit(thr', msDt(:, 3), f);
catch emsg
	fprintf(1, 'DN_OutlierInclude: error fitting linear trend to valid proportion: %s\n', emsg);
end

if isempty(emsg)
	out.nfla = c.a;
	out.nflb = c.b;
	out.nflr2 = gof.rsquare;
	out.nflrmse = gof.rmse;
else
	out.nfla = NaN;
	out.nflb = NaN;
	out.nflr2 = NaN;
	out.nflrmse = NaN;
end

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
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [5, 1, 15]);
f = fittype('a*exp(b*x)+c', 'options', s);
emsg = [];
try
	[c, gof] = fit(thr', msDt(:, 6), f);
catch emsg
	warning('Error fitting exponential growth to std: %s\n', emsg.message);
end

if isempty(emsg)
	out.stdrfexpa = c.a;
	out.stdrfexpb = c.b;
	out.stdrfexpc = c.c;
	out.stdrfexpr2 = gof.rsquare;
	out.stdrfexprmse = gof.rmse;
else
	out.stdrfexpa = NaN;
	out.stdrfexpb = NaN;
	out.stdrfexpc = NaN;
	out.stdrfexpr2 = NaN;
	out.stdrfexprmse = NaN;
end

% ------------------------------------------------------------------------------
%% Fit linear to errors in range
% ------------------------------------------------------------------------------
s = fitoptions('Method', 'NonlinearLeastSquares', 'StartPoint', [40, 4]);
f = fittype('a*x +b', 'options', s);
emsg = '';
try
	[c, gof] = fit(thr', msDt(:, 6), f);
catch emsg
	fprintf(1, 'DN_OutlierInclude: error fitting linear trend to std: %s\n', emsg);
end

if isempty(emsg)
	out.stdrfla = c.a;
	out.stdrflb = c.b;
	out.stdrflr2 = gof.rsquare;
	out.stdrflrmse = gof.rmse;
else
	out.stdrfla = NaN;
	out.stdrflb = NaN;
	out.stdrflr2 = NaN;
	out.stdrflrmse = NaN;
end

if doPlot
	figure('color', 'w')
	errorbar(thr, msDt(:, 1), msDt(:, 2), 'k');
end

end
