function out = MF_ExpSmoothing(x, ntrain, alpha)
% MF_ExpSmoothing   Exponential smoothing as a one-step forecaster: the best smoothing parameter and its residuals.
%
% Fits an exponential smoothing model to the time series, in which the forecast is an
% exponentially weighted average of past values, S(t) = alpha*X(t) + (1-alpha)*S(t-1),
% and alpha is the smoothing parameter. The best alpha is found by minimizing the
% root-mean-square one-step prediction error over a training set (the first ntrain
% samples): a coarse search over five values from 0.1 to 0.9, with a parabola fitted
% to the three lowest errors, and then a finer search around its minimum. The chosen
% alpha is then used to forecast each sample after the training set from the samples
% before it, and the outputs report alpha and statistics of the residuals of these
% held-out forecasts, from the shared residual summary MF_ResidualAnalysis (at its
% 'full' level).
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
% alphamin_1, the first estimate of it: the minimum of the parabola fitted in the
%           coarse search (not bounded to the interval 0 to 1)
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
	%% (*) Optimize alpha (*)
	% optimize alpha over the training set xtrain. This is a choice to only use
	% the first section of the time series. The length of the training set is
	% set in the input
	% optimize alpha based on the defined training set of time series values

	% Training set xtrain
	xtrain = x(1:ntrain);

	% Do a descent: fit a quadratic to available points.
	% (This replaces a dead 'bruteForce' alternative that was permanently disabled by a
	%  hardcoded flag and could not have run as written -- it referenced an undefined
	%  variable 'j' in 'rmses(j)', an undefined 'rmse_n', and called the interactive
	%  'input()', which would have hung a batch job. This descent branch was always the
	%  only code path actually exercised.)

	% (1) use alpha = 0.01, 0.1, 0.5, 0.8
	%         alphar = [0.1, 0.2, 0.5, 0.8, 0.9];
	alphar = linspace(0.1, 0.9, 5);
	rmses = zeros(4, 1);

	for k = 1:length(alphar)
		a = alphar(k);

		xf = SUB_fit_exp_smooth(xtrain, a);

		% Issue forecasts
		fore = xf(3:end);
		orig = xtrain(3:end);
		rmses(k) = sqrt(mean((fore - orig).^2)); % compute rmse
	end

	% fit quadratic to set alpha
	[sort_rmses, ix] = sort(rmses);
	rkeep = ix(1:3); % fit on 3 points closest to minimum
	p = polyfit(alphar(rkeep)', rmses(rkeep), 2);
	aar = (0:0.005:1);
	y = polyval(p, aar);
	%         plot(aar,y,':k'); hold on;
	%         plot(alphar,rmses,'or');
	%         plot(alphar(rkeep),rmses(rkeep),'*m'); hold off
	alphamin = -p(2) / (2 * p(1));
	out.alphamin_1 = alphamin;
	out.p1_1 = abs(p(1)); % concavity
	out.cup_1 = sign(p(1));

	if p(1) < 0 % concave down -- it's looking at a maximum
		% weird case
		if y(1) < y(end);
			alphamin = 0.01;
		else
			alphamin = 1;
		end
	else
		% Search again around this
		alphar = linspace(alphamin - 0.1, alphamin + 0.1, 5);
		if any(alphar <= 0)
			alphar = linspace(0.01, max(alphamin, 0) + 0.1, 5);
		elseif any(alphar >= 1)
			alphar = linspace(min(alphamin, 1) - 0.1, 1, 5);
		end

		for k = 1:length(alphar)
			a = alphar(k);

			xf = SUB_fit_exp_smooth(xtrain, a);

			% Issue forecasts
			fore = xf(3:end);
			orig = xtrain(3:end);
			rmses(k) = sqrt(mean((fore - orig).^2)); % compute rmse

		end

		% Fit quadratic to set alpha
		p = polyfit(alphar', rmses, 2);
		%             aar = 0:0.005:1;
		%             y = polyval(p,aar);
		%             plot(aar,y,':k'); hold on;
		%             plot(alphar,rmses,'or'); hold off

		if p(1) < 0
			alphamin = alphar(rmses == min(rmses));
			% This is quite bad -- the first step didn't find a local
			% minimum...
		else % minimum of quadratic fit
			alphamin = -p(2) / (2 * p(1));
			if alphamin > 1, alphamin = 1; end
			if alphamin <= 0, alphamin = 0.01; end
		end

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
	% The forecast of x(ii+1) restarts the smoother at the start of the series, with
	% initial value s(1) = mean(x(1:ii-1)), and runs s(jj) = a*x(jj) + (1-a)*s(jj-1)
	% for jj = 2:ii. Unrolling the recursion, the forecast is
	%   s(ii) = (1-a)^(ii-1)*mean(x(1:ii-1)) + sum_{jj=2}^{ii} a*(1-a)^(ii-jj)*x(jj),
	% where the sum is itself a one-pole filter of x(2:end), so all forecasts are
	% found in O(N) rather than by restarting the loop at every ii.
	x = x(:);
	nx = length(x);
	xf = zeros(nx, 1);
	if nx < 3
		return
	end

	ii = (2:nx - 1)';
	runMean = cumsum(x(1:nx - 2)) ./ (1:nx - 2)'; % mean(x(1:ii-1))
	ewma = filter(a, [1, -(1 - a)], [0; x(2:end)]); % sum_{jj=2}^{ii} a*(1-a)^(ii-jj)*x(jj)

	% S(t) = Xf(t) is forecasted value for X(t+1)
	xf(ii + 1) = (1 - a).^(ii - 1) .* runMean + ewma(ii);
end

end
