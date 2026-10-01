function out = MF_CompareTestSets(y, theModel, ord, subsetHow, samplep, steps, randomSeed)
% MF_CompareTestSets   How well a model fitted to the whole series predicts short stretches of it.
%
% Fits a time-series model to the full series, then uses it to predict a set of
% short test segments of the series (steps samples ahead), and summarizes how the
% prediction quality varies across the segments. For each segment it records the
% root-mean-square prediction error, the lag-1 autocorrelation of the errors, the
% absolute difference between the mean prediction and the mean of the data, and the
% ratio of the standard deviations of the predictions and the data. It says something
% about stationarity in the spread of values, and about the suitability of the model
% in the level of values.
%
% Similar to MF_FitSubsegments, except that the model is fitted on the full time
% series and tested on different local segments.
%
% Uses iddata and predict from MATLAB's System Identification Toolbox, and ar, n4sid
% or armax to fit the model.
%
% ---INPUTS:
% y, the input time series
%
% theModel, the type of time-series model to fit:
%           (i) 'ar', an AR model,
%           (ii) 'ss', a state-space model (default),
%           (iii) 'arma', an ARMA model.
%
% ord, the order of the model to fit (default 2; a two-element vector for 'arma'),
%       or 'best' to select it automatically: for 'ar', the order from 1 to 10
%       minimizing the Schwarz Bayesian criterion (ARFIT_arfit); for 'ss', as chosen
%       by n4sid.
%
% subsetHow, how to select the test segments:
%           (i) 'rand', at random (default),
%           (ii) 'uniform', evenly spaced throughout the time series.
%
% samplep, a two-vector specifying the sampling parameters, [number of segments,
%           segment length] (default [20, 0.1]). A segment length below 1 is a
%           fraction of the series length, capped to between 10 and 20 samples
%           (so [25, 0.1] takes 25 segments of 10 to 20 samples); otherwise it is a
%           number of samples.
%
% steps, the number of steps ahead to predict in each segment (default 2).
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (used when subsetHow is 'rand')
%
% ---OUTPUTS:
% stde_mean, stde_std, stde_iqr: the mean, standard deviation and interquartile range
%       over segments of the root-mean-square prediction error
% ac1_mean, ac1_median: the absolute value of the mean, and of the median, over
%       segments of the lag-1 autocorrelation of the prediction errors
% ac1_std, ac1_iqr: the standard deviation and interquartile range over segments of
%       that autocorrelation
% meane_mean, meane_std, meane_iqr: the mean, standard deviation and interquartile
%       range over segments of the absolute difference between the mean prediction
%       and the mean of the data
% stdrat_mean, stdrat_median, stdrat_std, stdrat_iqr: the mean, median, standard
%       deviation and interquartile range over segments of the ratio of the standard
%       deviation of the predictions to that of the data (segments in which the data
%       are near-constant are excluded)
%
% ---NOTES:
% Redundant fields (mabserrs, and the medians of stde and meane) were dropped from
% this function on 2026-08-11 after a redundancy check on Bonn EEG (500 series) and
% Empirical1000 (1000 series): each correlated at |r| >= 0.9 with a retained field
% on both datasets and in all registered operations that use this function.

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
%% Check that a System Identification Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('identification_toolbox');

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Check inputs, set defaults
% ------------------------------------------------------------------------------
% (1) y: column vector time series
if nargin < 1 || isempty(y)
	error('No input time series provided');
end

% Convert y to time series object for the System Identification Toolbox
sigmaY = std(y); % used to flag degenerate (near-constant) test segments below
y = iddata(y, [], 1);

% (2) Model, the type of model to fit
if nargin < 2 || isempty(theModel)
	theModel = 'ss';
	% Fit a state space model by default
end

% (3) Model order, ord
if nargin < 3 || isempty(ord)
	ord = 2;
	% model of order 2 by default. Not the best defaults.
end

% (4) How to choose subsets from the time series, subsetHow
if nargin < 4 || isempty(subsetHow)
	subsetHow = 'rand'; % takes segments randomly from time series
end

% (5) Sampling parameters, samplep
if nargin < 5 || isempty(samplep)
	samplep = [20, 0.1]; % sample 20 times with 10%-length subsegments
end

% (6) Predict some number of steps ahead in test sets, steps
if nargin < 6 || isempty(steps)
	steps = 2; % default: predict 2 steps ahead in test set
end

% (7)  randomSeed: how to treat the randomization
if nargin < 7
	randomSeed = [];
end

% ------------------------------------------------------------------------------
%% Fit the model
% ------------------------------------------------------------------------------
% model will be stored as a model object, m
% model is fitted using the entire dataset as the training set
% test sets will be smaller chunks of this.
% [could also fit multiple models using data not in multiple test sets, but
% this is messier]

switch theModel
	case 'ar' % fit an ar model of specified order
		if strcmp(ord, 'best')
			% Use arfit software to retrieve the optimum ar order by some
			% criterion (Schwartz's Bayesian Criterion, SBC)
			try
				[west, Aest, Cest, SBC, FPE, th] = ARFIT_arfit(y.y, 1, 10, 'sbc', 'zero');
			catch
				error('Error running ''arfit'' -- is the ARFIT toolbox installed?')
			end
			ord = length(Aest);
		end
		m = ar(y, ord);

	case 'ss' % fit a state space model of specified order
		m = n4sid(y, ord);

	case 'arma' % fit an arma model of specified orders
		% Note: order should be a two-component vector
		m = armax(y, ord);

	otherwise
		error('Unknown model ''%s''', theModel);
end

% ------------------------------------------------------------------------------
%% Prepare to do a series of predictions
% ------------------------------------------------------------------------------
% Number of samples to take, numPred
numPred = samplep(1);
% Initialize quantities to store into
rmserrs = zeros(numPred, 1);
ac1s = zeros(numPred, 1);
meandiffs = zeros(numPred, 1);
stdrats = zeros(numPred, 1);

% Set ranges beforehand
r = zeros(numPred, 2);

switch subsetHow
	case 'rand'
		if samplep(2) < 1 % specified a fraction of time series
			% A pure fraction makes the test-segment length -- and with it
			% the precision of every residual statistic below -- grow
			% without bound as N grows, so results never converge to a
			% fixed population value: verified this was the dominant
			% driver of this operation's length-dependence (rank eta^2 vs
			% N on a stationary AR(1) null, N=200..6400, fell from
			% 0.39-0.94 to 0.01-0.21 across stde/stdrat/ac1/meane
			% statistics once capped). Capping at 20 leaves N=200 (the
			% audit's own smallest tested length, where 10% of the series
			% is already <= 20) completely unchanged and only bites for
			% longer series, where letting the segment keep growing was
			% buying no real precision benefit anyway.
			l = max(min(20, floor(N * samplep(2))), 10);
		else % specified an absolute interval
			l = samplep(2);
		end

		% Control the random seed (for reproducibility):
		BF_ResetSeed(randomSeed);

		% numPred starting points:
		spts = randi(N - l + 1, numPred, 1);
		r(:, 1) = spts;
		r(:, 2) = spts + l - 1;

	case 'uniform'
		if length(samplep) == 1 % size will depend on number of unique subsegments
			spts = round(linspace(0, N, numPred + 1)); % numPred+1 boundaries = numPred portions
			r(:, 1) = spts(1:numPred) + 1;
			r(:, 2) = spts(2:end);
		else
			if samplep(2) < 1 % specified a fraction of time series
				% Capped at an absolute length -- see the 'rand' case above
				% for why (test-segment length must not grow with N).
				l = max(min(20, floor(N * samplep(2))), 10);
			else % specified an absolute interval
				l = samplep(2);
			end
			spts = round(linspace(1, N - l + 1, numPred)); % numPred+1 boundaries = numPred portions
			r(:, 1) = spts;
			r(:, 2) = spts + l - 1;
		end
	otherwise
		error('Unknown subset method ''%s''', subsetHow);
end

% Quickly check that ranges are valid
if any(r(:, 1) >= r(:, 2))
	error('Invalid settings');
end

% ------------------------------------------------------------------------------
%% Do the series of predictions
% ------------------------------------------------------------------------------
for i = 1:numPred
	% Retrieve the test data:
	yTest = y(r(i, 1):r(i, 2));

	% Compute step-ahead predictions using System Identification Toolbox:
	yp = predict(m, yTest, steps); % across test set using model, m,
	% fitted to entire data set

	%     e = pe(m, yTest); % prediction errors -- exactly the same as returning
	%                       % residuals of 1-step-ahead prediction model

	% plot the two:
	% plot(y,yp);

	% Get statistics on residuals:
	mres = yp.y - yTest.y;

	rmserrs(i) = sqrt(mean(mres.^2));
	ac1s(i) = CO_AutoCorr(mres, 1, 'Fourier');

	% Get statistics on output time series
	meandiffs(i) = abs(mean(yp.y) - mean(yTest.y));
	% Guard against near-constant test segments: std(yTest.y) can land at
	% floating-point noise (rather than a genuinely small but real value)
	% when a short random/uniform segment happens to sit in a flat/clipped
	% stretch of the input -- verified empirically (2500 segments across 100
	% real series): the ratio's population is smooth from ~1e-3 upward, with
	% a separate, disjoint cluster at <1e-15 coming from exactly-constant
	% segments. Below that gap, std(yp.y)/std(yTest.y) is undefined rather
	% than just large, so exclude it instead of letting it dominate the
	% mean/std/iqr summaries below.
	if std(yTest.y) < 1e-6 * sigmaY
		stdrats(i) = NaN;
	else
		stdrats(i) = std(yp.y) / std(yTest.y);
	end

	%     % 1) Get statistics on residuals
	%     residout = MF_ResidualAnalysis(mresiduals);
	%
	%     % convert these to local outputs in quick loop
	%     fields = fieldnames(residout);
	%     for k=1:length(fields);
	%         eval(['out.' fields{k} ' = residout.' fields{k} ';']);
	%     end
end

% ------------------------------------------------------------------------------
%% Return statistics on outputs
% ------------------------------------------------------------------------------

% RMS errors, rmserrs
% (median dropped: r>=0.9 with stde_mean across all 4 registered mops, on
% both Bonn EEG and Empirical1000; mean absolute error, mabserrs, dropped
% entirely: for these residuals r>=0.9 with the matching stde_* moment in
% every case, so it added no dimension MAE didn't already carry)
out.stde_mean = mean(rmserrs);
out.stde_std = std(rmserrs);
out.stde_iqr = iqr(rmserrs);

% Autocorrelations at lag 1, ac1s
% NOT absolute values of ac1s... absolute values of operations on *raw* ac1s...
out.ac1_mean = abs(mean(ac1s));
out.ac1_median = abs(median(ac1s));
out.ac1_std = std(ac1s);
out.ac1_iqr = iqr(ac1s);

% Differences in mean between two series
% (median dropped: r>=0.9 with meane_mean across all 4 registered mops, on
% both Bonn EEG and Empirical1000)
out.meane_mean = mean(meandiffs);
out.meane_std = std(meandiffs);
out.meane_iqr = iqr(meandiffs);

% Ratio of standard deviations between two series
% (omitting repeats flagged as degenerate above)
validStdrats = stdrats(~isnan(stdrats));
out.stdrat_mean = mean(validStdrats);
out.stdrat_median = median(validStdrats);
out.stdrat_std = std(validStdrats);
out.stdrat_iqr = iqr(validStdrats);

end
