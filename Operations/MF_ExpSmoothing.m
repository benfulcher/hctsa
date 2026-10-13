function out = MF_ExpSmoothing(x, ntrain, alpha)
% MF_ExpSmoothing   Exponential smoothing as a one-step forecaster: the best smoothing parameter and its residuals.
%
% Fits a simple exponential smoothing model to the time series, in which the forecast
% of X(t) is the level S(t-1), an exponentially weighted average of past values,
% S(t) = alpha*X(t) + (1-alpha)*S(t-1), with the initial level the mean of the first
% (up to) 10 samples; alpha is the smoothing parameter. The best alpha is the value in
% [0.01, 1] that minimizes the root-mean-square one-step prediction error over a
% training set (the first ntrain samples, all but the first forecast scored): the error
% is evaluated on the grid 0.01, 0.02, ..., 1 and the best grid point g is refined by 40
% golden-section steps on [max(0.01, g - 0.01), min(1, g + 0.01)]. A coarse search over
% five values from 0.1 to 0.9, with a parabola fitted to the three lowest errors, is
% also reported. The chosen alpha is then used to forecast each sample after the
% training set from the samples before it, and the outputs report alpha and statistics
% of the residuals of these held-out forecasts, from the shared residual summary
% MF_ResidualAnalysis (at its 'full' level).
%
% ---INPUTS:
% x, the input time series
%
% ntrain, the number of samples to use for training (can be a proportion of the
%           time-series length, if between 0 and 1). It is kept between 100 and 1000
%           samples. Default is min(100, N). If the series is shorter than ntrain, or
%           fewer than 50 samples remain after the training set, a NaN is returned.
%
% alpha, the exponential smoothing parameter, or 'best' (default) to fit it on the
%           training set.
%
% ---OUTPUTS (when alpha is 'best'; the residual statistics are those of the
% held-out samples, after the first ntrain):
% alphamin, the fitted smoothing parameter (between 0.01 and 1)
% slope05_1, the slope at alpha = 0.5 of the parabola fitted to the three lowest errors
%           of the coarse search (negative if the error falls towards larger alpha)
% p1_1, the size of the quadratic coefficient of that parabola
% cup_1, the sign of that coefficient (+1 if the parabola opens upward)
% meane, mean of the residuals (prediction minus data)
% meanabs, mean absolute residual
% stde, standard deviation of the residuals
% maxonstd, largest absolute residual, in units of the residual standard deviation
% ac1, ac2, ac3: autocorrelation of the (z-scored) residuals at lags 1, 2 and 3
% propbth, proportion of the residual autocorrelations at lags 1 to 25 within the
%           significance band +/- 2.6/sqrt(N)
% taurat, decorrelation time of the residuals (first zero-crossing of their
%           autocorrelation function) divided by that of the time series
% ftbth, first lag at which the residual autocorrelation falls inside the
%           significance band (26 if it never does)
% normksstat, Kolmogorov-Smirnov statistic of the residuals against a Gaussian
% sws, standard deviation across 5 windows of the local standard deviation of the
%           residuals, relative to their overall standard deviation
% swm, standard deviation across 5 windows of the local mean of the residuals,
%           relative to their overall standard deviation
% popt, the order (1 to 10) of the AR model fitted to the residuals, selected by the
%           Schwarz Bayesian criterion
% minsbc, the corresponding Schwarz Bayesian criterion
%
% ---REFERENCES:
% C. Chatfield, "The Analysis of Time Series", CRC Press LLC (2004).
% R. G. Brown, "Smoothing, Forecasting and Prediction of Discrete Time Series",
% Prentice-Hall (1963).
% J. Kiefer, "Sequential minimax search for a maximum", Proc. Amer. Math. Soc. 4, 502
% (1953) (golden-section search).
%
% ---NOTES:
% The residuals are those of the one-step forecasts of the held-out samples only
% (samples ntrain+1 to N), so they are not biased by the fit of alpha to the
% training set. Each forecast uses the samples before it, including earlier
% held-out ones. For the registered call (ntrain = 0.5) at N = 1000, 500 samples
% are held out.
%
% Code is adapted from that provided by Siddharth Arora (Siddharth.Arora@sbs.ox.ac.uk).

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

doPlot = false; % plot outputs
N = length(x); % the length of the time series

% -------------------------------------------------------------------------------
% Check Inputs:
% -------------------------------------------------------------------------------

% (*) ntrain -- either the number or proportion of training points
if nargin < 2 || isempty(ntrain)
	ntrain = min(100, N); % if the time series is shorter than 100 samples(!)
end
% Can give training set length as a proportion of the time-series length:
if (ntrain > 0) && (ntrain < 1)
	ntrain = floor(N * ntrain);
end

% Check training set size is between the following range:
minTrain = 100;
maxTrain = 1000;

if ntrain > maxTrain; % larger than maximum training set size
	fprintf(1, 'Training set size reduced from %u to maximum of 1000 samples.\n', ntrain);
	ntrain = 1000;
end
if ntrain < minTrain; % smaller than minimum training set size
	fprintf(1, 'Training set size increased from %u to minimum of 100.\n', ntrain);
	ntrain = 100;
end

if N < ntrain + 50 % too few samples held out after the training set
	fprintf(1, 'Time Series too short for exponential smoothing\n')
	out = NaN; return
end

% (*) alpha, the smoothing parameter
if nargin < 3 || isempty(alpha)
	alpha = 'best';
end

% -------------------------------------------------------------------------------
% Exponential smoothing
% Using: S(t+1) = A.X(t+1) + (1-A).S(t) where S(1) = X(1)
% Finding optimum parameter A [0,1] using RMSE

if strcmp(alpha, 'best')
	xtrain = x(1:ntrain);
	% Coarse search, reported through the parabola fitted to its three lowest errors:
	alphar = linspace(0.1, 0.9, 5);
	rmses = zeros(5, 1);
	for k = 1:length(alphar)
		rmses(k) = SUB_rmse(xtrain, alphar(k));
	end
	if any(~isfinite(rmses))
		out = NaN; return
	end
	[~, ix] = sort(rmses);
	rkeep = ix(1:3); % fit on the 3 points closest to the minimum
	p = polyfit(alphar(rkeep)', rmses(rkeep), 2);
	% Slope of the parabola at the middle of the scan, a + b for a*alpha^2 + b*alpha + c.
	% (Its vertex, -b/(2a), was reported as alphamin_1 before: a ratio that is unbounded
	% when the parabola is nearly flat. Together with p1_1 and cup_1 the slope carries the
	% same information, vertex = 0.5 - slope/(2a), and stays finite and continuous.)
	out.slope05_1 = p(1) + p(2);
	out.p1_1 = abs(p(1)); % concavity
	out.cup_1 = sign(p(1));

	% The best alpha: a grid of step 0.01, then golden-section refinement
	grid = (1:100)' / 100;
	gridErr = zeros(100, 1);
	for k = 1:100
		gridErr(k) = SUB_rmse(xtrain, grid(k));
	end
	[bestErr, ix] = min(gridErr);
	g = grid(ix);
	lo = max(0.01, g - 0.01); hi = min(1, g + 0.01);
	ratio = (sqrt(5) - 1) / 2;
	c = hi - ratio * (hi - lo); d = lo + ratio * (hi - lo);
	fc = SUB_rmse(xtrain, c); fd = SUB_rmse(xtrain, d);
	for k = 1:40
		if fc <= fd
			hi = d; d = c; fd = fc;
			c = hi - ratio * (hi - lo);
			fc = SUB_rmse(xtrain, c);
		else
			lo = c; c = d; fc = fd;
			d = lo + ratio * (hi - lo);
			fd = SUB_rmse(xtrain, d);
		end
	end
	aRef = (lo + hi) / 2;
	if SUB_rmse(xtrain, aRef) < bestErr
		alphamin = aRef;
	else
		alphamin = g;
	end
	out.alphamin = alphamin;
	alpha = alphamin;
end

if isnan(alpha)
	error('Alpha is a NaN?!')
end

%% (2) Fit to the whole time series

% % Plot in-sample error as a function of smoothing parameter
% figure(1);
% plot(alphar,rmses,'*');
% xlabel('\alpha - Smoothing parameter');
% ylabel('RMSE');

% Plot original time series and smoothed data using optimum values
y = SUB_fit_exp_smooth(x, alpha);

% Residuals only on the held-out part (after the ntrain samples used to fit alpha):
yp = y(ntrain+1:N); % predicted
xp = x(ntrain+1:N); % original
e = yp - xp; % residuals
% in_sample_error = sqrt(mean((yp-xp).^2));
% out.insamplermse = in_sample_error;

% -------------------------------------------------------------------------------
% Get statistics on residuals using MF_ResidualAnalysis
residout = MF_ResidualAnalysis(e, xp, 'full');

% Convert these to local outputs in quick loop:
fields = fieldnames(residout);
for k = 1:length(fields)
	out.(fields{k}) = residout.(fields{k});
end

if doPlot
	figure('color', 'w'); box('on')
	t = 1:length(yp);
	plot(t, xp, 'b', t, yp, 'k');
	legend('Obs', 'Fit');
	xlabel('Time');
	ylabel('Amplitude');
end

% ------------------------------------------------------------------------------
function xf = SUB_fit_exp_smooth(x, a)
	% One-step forecasts: xf(t) = S(t-1), with S(0) the mean of the first (up to) 10 samples
	x = x(:);
	level0 = mean(x(1:min(10, length(x))));
	level = filter(a, [1, -(1 - a)], x, (1 - a) * level0);
	xf = [level0; level(1:end-1)];
end
% ------------------------------------------------------------------------------
function e = SUB_rmse(x, a)
	% Root-mean-square one-step error, all but the first forecast
	xf = SUB_fit_exp_smooth(x, a);
	e = sqrt(mean((xf(2:end) - x(2:end)).^2));
end

end
