function out = SC_FluctAnal(x, q, wtf, tauStep, k, lag, logInc)
% SC_FluctAnal   Scaling exponents from fluctuation analysis, by a choice of methods.
%
% The series is integrated (cumulative sum) and cut into windows of tau samples,
% a measure of the size of the fluctuations is taken in each window, and the
% fluctuation function F(tau) is the q-th order mean of these (q = 2 is the
% usual root-mean-square). For a self-similar series F(tau) ~ tau^alpha, so the
% scaling exponent alpha is the slope of log F(tau) against log tau. The methods
% differ in how fluctuations are measured in a window (wtf, below). Much of the
% implementation follows the discussion of scaling methods in Talkner and Weber
% (2000).
%
% Timescales tau run from 5 samples to half the series length, in tauStep
% logarithmically spaced steps (logInc true) or in steps of tauStep samples
% (logInc false, in which case log F is interpolated by a spline onto 50 evenly
% spaced values of log tau). Peng et al. (1995) used 5 to a quarter of the series
% length; Little et al. (2007) used 4 to half. Fewer than 8 timescales give NaN.
%
% ---INPUTS:
% x, the input time series
% q, the order of the fluctuation function (2, the usual choice, gives
%       root-mean-square fluctuations); default 2
% wtf, what to fluctuate (default 'rsrange'):
%       'dfa': subtract a polynomial trend of order k in each window, then take
%           the q-th order mean over all points (detrended fluctuation analysis)
%       'endptdiff': the difference between the end points of each window
%       'range': the range in each window
%       'std': the standard deviation in each window (cf. Cannon et al. 1997)
%       'iqr': the interquartile range in each window
%       'rsrange': the range after subtracting the straight line joining the end
%           points of the window (cf. Caccia et al. 1997)
%       'rsrangefit': the range after subtracting a polynomial trend of order k
%       'nothing': no windowing statistic; the q-th order mean over all points of
%           the integrated series that fit in whole windows, with no detrending
% tauStep, the number of timescales (logInc true) or the step in tau, in samples
%       (logInc false); default 1
% k, the polynomial order of the detrending, for 'dfa' and 'rsrangefit'; default 1
% lag, an optional time lag for the integrated profile (Alvarez-Ramirez et al.
%       2009): if given, the series is subsampled by taking every lag-th point
%       before integrating; default none
% logInc, whether the timescales are logarithmically spaced (true, the
%       recommended choice) or linearly spaced (false); default true
%
% ---OUTPUTS: statistics of a robust linear fit of log F(tau) against log tau, and
% of fitting two straight lines to the same data, with the split point chosen to
% minimize the combined fitting error (each line spans at least a quarter of the
% timescales, and at least 8):
% linfitint, alpha, se1, se2, ssr, resac1: the intercept, slope (the scaling
%       exponent alpha), standard errors of the intercept and the slope, mean
%       squared residual, and lag-1 autocorrelation of the residuals, of the
%       single line fitted over all timescales
% r1_linfitint, r1_alpha, r1_se1, r1_se2, r1_ssr, r1_resac1: the same, for the
%       first (shorter-timescale) line of the two-line fit
% r2_linfitint, r2_alpha, r2_se1, r2_se2, r2_ssr, r2_resac1: the same, for the
%       second (longer-timescale) line
% logtausplit, the value of log(tau) at the split between the two lines
% prop_r1, the proportion of the timescales covered by the first line
% ratsplitminerr, the ratio of the minimum two-line fitting error (mean squared
%       error pooled over both lines) to ssr
% meanssr, stdssr, the mean and the standard deviation of the two-line fitting
%       error across the candidate split points
% alpharat, the ratio r1_alpha / r2_alpha
% (All are NaN if there are too few timescales; the two-line fields are NaN if
% the timescales are too few to support two lines.)
%
% ---REFERENCES:
% P. Talkner and R. O. Weber, "Power spectrum and detrended fluctuation analysis:
% Application to daily temperatures", Phys. Rev. E 62(1), 150 (2000).
% M. J. Cannon et al., "Evaluating scaled windowed variance methods for estimating
% the Hurst coefficient of time series", Physica A 241(3-4), 606 (1997).
% D. C. Caccia et al., "Analyzing exact fractal time series: evaluating
% dispersional analysis and rescaled range methods", Physica A 246(3-4), 609 (1997).
% J. Alvarez-Ramirez et al., "Using detrended fluctuation analysis for lagged
% correlation analysis of nonstationary signals", Phys. Rev. E 79(5), 057202 (2009).
% C.-K. Peng et al., "Statistical properties of DNA sequences", Physica A
% 221(1-3), 180 (1995).
% Little et al., "Exploiting Nonlinear Recurrence and Fractal Scaling Properties
% for Voice Disorder Detection", Biomed. Eng. Online 6, 23 (2007).
%
% ---NOTES:
% In hctsa the function is also applied to the absolute values of the z-scored
% series and to its sign (+1 above the mean, -1 below), to look at the scaling of
% the magnitude and of the sign of the fluctuations separately.

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
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(q)
	q = 2; % RMS fluctuations
end
if nargin < 3 || isempty(wtf)
	wtf = 'rsrange'; % re-scaled range analysis by default
end
if nargin < 4 || isempty(tauStep)
	% the increment of tau (for linear)
	% or number of points in logarithmic range (for logarithmic)
	tauStep = 1;
end
if nargin < 5 || isempty(k)
	k = 1; % often not needed, only for 'dfa' and 'rsrangefit'
end
if nargin < 6
	lag = '';
end
if nargin < 7
	logInc = true;
end

% ------------------------------------------------------------------------------
N = length(x); % length of the time series
doPlot = false; % plot relevant outputs to figure

% -------------------------------------------------------------------------------
% 1) Compute integrated sequence
if isempty(lag) || lag == 1
	% A normal cumsum:
	y = cumsum(x);
else
	% If a lag is specified, do a decimation:
	y = cumsum(x(1:lag:end));
end

% -------------------------------------------------------------------------------
% Perform scaling over a range of tau, up to a fifth the time-series length
% -------------------------------------------------------------------------------
% Peng (1995) suggests 5:N/4 for DFA
% Caccia suggested from 10 to (N-1)/2...
% -------------------------------------------------------------------------------
if logInc
	taur = unique(round(exp(linspace(log(5), log(floor(N / 2)), tauStep))));
	% in this case tauStep is the number of points to compute
else
	taur = 5:tauStep:floor(N / 2); % maybe increased??
end
ntau = length(taur); % analyze the time series across this many timescales

if ntau < 8 % fewer than 8 points
	fprintf(1, 'This time series (N = %u) is too short to analyze using this fluctuation analysis\n', N);
	out = NaN;
	return
end

% -------------------------------------------------------------------------------
% 2) Compute the fluctuation function, F
F = zeros(1, ntau);
% Each entry correponds to a given scale, tau

for i = 1:ntau
	% buffer the time series at the scale tau
	tau = taur(i); % the scale on which to compute fluctuations

	y_buff = buffer(y, tau);
	if size(y_buff, 2) > floor(N / tau) % zero-padded, remove trailing set of points...
		y_buff = y_buff(:, 1:end - 1);
	end

	% analyzed length of time series (with trailing end-points removed)
	nn = size(y_buff, 2) * tau;

	switch wtf
		case 'nothing'
			y_dt = reshape(y_buff, nn, 1);
		case 'endptdiff'
			% look at differences in end-points in each subsegment
			y_dt = y_buff(end, :) - y_buff(1, :);
		case 'range'
			y_dt = range(y_buff);
		case 'std'
			% something like what they do in Cannon et al., Physica A 1997,
			% except with bridge/linear detrending and overlapping segments
			% (scaled windowed variance methods). But I think we have
			% enough of this sort of thing already...
			y_dt = std(y_buff);
		case 'iqr'
			y_dt = iqr(y_buff);
		case 'dfa'
			tt = (1:tau)'; % faux time range
			for j = 1:size(y_buff, 2);
				% fit a polynomial of order k in each subsegment
				p = polyfit(tt, y_buff(:, j), k);
				% remove the trend, store back in y_buff
				y_buff(:, j) = y_buff(:, j) - polyval(p, tt);
			end
			% reshape to a column vector, y_dt (detrended)
			y_dt = reshape(y_buff, nn, 1);
		case 'rsrange'
			% Remove straight line first: Caccia et al. Physica A, 1997
			% Straight line connects end points of each window:
			b = y_buff(1, :);
			m = y_buff(end, :) - b;
			y_buff = y_buff - (linspace(0, 1, tau)' * m + ones(tau, 1) * b);
			y_dt = range(y_buff);
		case 'rsrangefit' % polynomial fit (order k) rather than endpoints fit: (~DFA)
			tt = (1:tau)'; % faux time range
			for j = 1:size(y_buff, 2);
				% fit a polynomial of order k in each subsegment
				p = polyfit(tt, y_buff(:, j), k);
				% remove the trend, store back in y_buff
				y_buff(:, j) = y_buff(:, j) - polyval(p, tt);
			end
			y_dt = range(y_buff);
		otherwise
			error('Unknown fluctuation analysis method ''%s''', wtf);
	end

	% Compute fluctuation function:
	F(i) = (mean(y_dt.^q)).^(1 / q);
end

% -------------------------------------------------------------------------------
% Smooth unevenly-distributed points in log space:
% -------------------------------------------------------------------------------
if logInc
	logtt = log(taur);
	logFF = log(F);
	numTimeScales = ntau;
else % need to smooth the unevenly-distributed points (using a spline)
	logtaur = log(taur); logF = log(F);
	numTimeScales = 50; % number of sampling points across the range
	logtt = linspace(min(logtaur), max(logtaur), numTimeScales); % even sampling in tau
	logFF = spline(logtaur, logF, logtt);
end

% -------------------------------------------------------------------------------
% Linear fit the log-log plot: full range
% -------------------------------------------------------------------------------
out = struct();
out = DoRobustLinearFit(out, logtt, logFF, 1:numTimeScales, '');

% PLOT THIS?:
if doPlot
	figure('color', 'w');
	plot(logtt, logFF, 'o-k');
	title(out.alpha)
end

%% WE NEED SOME SORT OF AUTOMATIC DETECTION OF GRADIENT CHANGES/NUMBER
%% OF PIECEWISE LINEAR PIECES

% ------------------------------------------------------------------------------
%% Try assuming two components (2 distinct scaling regimes)
% ------------------------------------------------------------------------------
% Move through, and fit a straight line to loglog before and after each point.
% Find point with the minimum sum of squared errors

% First spline interpolate to get an even sampling of the interval
% (currently, in the log scale, there are relatively more at slower timescales)

% Determine the errors
% ------------------------------------------------------------------------------
% minPoints scales with numTimeScales rather than being a fixed constant: a
% fixed small minPoints (e.g. 6) lets the search reach breakpoints right at the
% edge of the domain, where a segment of only a handful of points can trivially
% achieve near-zero fit error regardless of whether the series has any real
% change in scaling behaviour. Confirmed empirically: for monofractal fGn (no
% true crossover), the raw fit-error curve is monotonic across the *entire*
% search domain, so the "best" split is always whichever end the search is
% allowed to reach, not a genuine interior minimum. Requiring each segment to
% span at least a quarter of the timescale range keeps the search away from
% these degenerate edge solutions. The floor of 8 also guarantees
% DoRobustLinearFit's own length>=8 requirement is always satisfied whenever a
% fit is attempted here, so the search bound and the fit's minimum-length
% requirement can no longer disagree (previously minPoints=6 could select a
% 6- or 7-point segment that DoRobustLinearFit would then discard as NaN).
sserr = nan(numTimeScales, 1); % don't choose the end points
minPoints = max(8, round(0.25 * numTimeScales));
if numTimeScales >= 2 * minPoints
	for i = minPoints:numTimeScales - minPoints
		r1 = 1:i;
		p1 = polyfit(logtt(r1), logFF(r1), 1);
		r2 = i:numTimeScales;
		p2 = polyfit(logtt(r2), logFF(r2), 1);
		% Mean squared error, pooled across both segments and normalized by
		% the total number of points sampled (numTimeScales), so that
		% ratsplitminerr below (which divides by out.ssr, a mean squared
		% error) is a genuinely comparable, tauStep-invariant ratio. This was
		% previously a straight sum of L2 norms (unnormalized), which scales
		% with sqrt(#points) and made ratsplitminerr roughly triple just from
		% varying tauStep 20->200 on an identical series -- not a real effect,
		% purely a units mismatch.
		e1 = polyval(p1, logtt(r1)) - logFF(r1);
		e2 = polyval(p2, logtt(r2)) - logFF(r2);
		sserr(i) = (sum(e1.^2) + sum(e2.^2)) / numTimeScales;
	end
end

if all(isnan(sserr))
	% Too few timescales to fit two distinct linear regimes meaningfully
	r1 = []; r2 = [];
	out.prop_r1 = NaN;
	out.logtausplit = NaN;
	out.ratsplitminerr = NaN;
	out.meanssr = NaN;
	out.stdssr = NaN;
else
	% breakPt is the point where it's best to fit a line before and another line after
	breakPt = find(sserr == min(sserr), 1, 'first');
	r1 = 1:breakPt;
	r2 = breakPt:numTimeScales;

	% Proportion of the domain of timescales corresponding to the first good linear fit
	out.prop_r1 = length(r1) / numTimeScales;

	out.logtausplit = logtt(breakPt);
	out.ratsplitminerr = min(sserr) / out.ssr;
	out.meanssr = nanmean(sserr);
	out.stdssr = nanstd(sserr);
end

if doPlot
	subplot(3, 1, 1)
	plot(y)
	subplot(3, 1, 2)
	plot(logtt(r1), logFF(r1), 'o-b')
	hold on;
	plot(logtt(r2), logFF(r2), 'o-r')
	subplot(3, 1, 3)
	plot(logtt, sserr, 'x-k')
end

% Check that at least 3 points are available

% -------------------------------------------------------------------------------
% Now we perform the robust linear fitting and get statistics on the two segments
% -------------------------------------------------------------------------------
% R1:
out = DoRobustLinearFit(out, logtt, logFF, r1, 'r1_');

% R2:
out = DoRobustLinearFit(out, logtt, logFF, r2, 'r2_');

if isnan(out.r1_alpha) || isnan(out.r2_alpha)
	out.alpharat = NaN;
else
	out.alpharat = out.r1_alpha / out.r2_alpha;
end

% -------------------------------------------------------------------------------
function out = DoRobustLinearFit(out, logtt, logFF, theRange, fieldName)
	% Get robust linear fit statistics on scaling range
	% Adds fields to the output structure

	if length(theRange) < 8 || all(isnan(logFF(theRange)))
		out.([fieldName, 'linfitint']) = NaN;
		out.([fieldName, 'alpha']) = NaN;
		out.([fieldName, 'se1']) = NaN;
		out.([fieldName, 'se2']) = NaN;
		out.([fieldName, 'ssr']) = NaN;
		out.([fieldName, 'resac1']) = NaN;
	else
		[linfit, stats] = robustfit(logtt(theRange), logFF(theRange));

		out.([fieldName, 'linfitint']) = linfit(1); % linear fit intercept
		out.([fieldName, 'alpha']) = linfit(2); % linear fit gradient
		out.([fieldName, 'se1']) = stats.se(1); % standard error in intercept
		out.([fieldName, 'se2']) = stats.se(2); % standard error in mean
		out.([fieldName, 'ssr']) = mean(stats.resid.^2); % mean squares residual
		out.([fieldName, 'resac1']) = CO_AutoCorr(stats.resid, 1, 'Fourier');
	end
end

end
