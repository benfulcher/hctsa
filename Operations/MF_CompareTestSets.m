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
% The model is fitted with BF_FitPredictorModel (AR: forward-backward least squares;
% state space: N4SID subspace identification; ARMA: conditional least squares) and its
% predictions of each segment start from a diffuse (completely unknown) initial state
% (BF_SSPredict).
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
%       minimizing the Schwarz Bayesian criterion (ARFIT_arfit); for 'ss', the order
%       from 1 to 10 minimizing Akaike's information criterion (BF_SubspaceID).
%
% subsetHow, how to select the test segments:
%           (i) 'rand', at random (default),
%           (ii) 'uniform', evenly spaced throughout the time series,
%           (iii) 'spread', start points spread over the series by the golden-ratio
%                 (Weyl) sequence: deterministic, evenly spread, and without the
%                 aliasing a regular spacing has with periodic series.
%
% samplep, a two-vector specifying the sampling parameters, [number of segments,
%           segment length] (default [20, 0.1]). A segment length below 1 is a
%           fraction of the series length, capped to between 10 and 20 samples
%           (so [25, 0.1] takes 25 segments of 10 to 20 samples); otherwise it is a
%           number of samples.
%
% steps, the number of steps ahead to predict in each segment (default 2).
%
% randomSeed, the seed of the random start points (see BF_RandomSeed; the numbers come
%               from the portable generator BF_Random; used when subsetHow is 'rand'). The
%               segments of a random draw are what the outputs depend on most, so
%               prefer 'spread' or 'uniform'.
%
% ---OUTPUTS:
% stde_mean, stde_std, stde_iqr: the mean, standard deviation and interquartile range
%       over segments of the root-mean-square prediction error
% ac1_mean, ac1_median: the absolute value of the mean, and of the median, over
%       segments of the lag-1 autocorrelation of the prediction errors
% ac1_std, ac1_iqr: the standard deviation and interquartile range over segments of
%       that autocorrelation (segments whose prediction errors do not vary beyond the
%       level of rounding, with a standard deviation at most 1e-9 times the larger of that
%       of the series and their mean absolute value, are excluded from the four ac1
%       statistics)
% meane_mean, meane_std, meane_iqr: the mean, standard deviation and interquartile
%       range over segments of the absolute difference between the mean prediction
%       and the mean of the data
% stdrat_mean, stdrat_median, stdrat_std, stdrat_iqr: the mean, median, standard
%       deviation and interquartile range over segments of the ratio of the standard
%       deviation of the predictions to that of the data (segments in which the data
%       are near-constant are excluded)
%
% ---NOTES:
% The first samples of each segment determine the initial state of the predictor and
% are not predicted (the first n + steps - 1 for a model of order n); the statistics of a
% segment use the samples predicted, and are NaN for a segment with fewer than 3 of them.
% All outputs are NaN if the model cannot be fitted.
% Redundant fields (mabserrs, and the medians of stde and meane) are not returned:
% each is largely redundant with a retained field.

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
N = length(y); % length of time series

% ------------------------------------------------------------------------------
%% Check inputs, set defaults
% ------------------------------------------------------------------------------
% (1) y: column vector time series
if nargin < 1 || isempty(y)
	error('No input time series provided');
end

sigmaY = std(y); % used to flag degenerate (near-constant) test segments below
y = y(:);

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

try
	[A, K, C] = BF_FitPredictorModel(y, theModel, ord);
catch emsg
	if contains(emsg.message, 'Unknown model')
		rethrow(emsg);
	end
	out = NaN; return % the model cannot be fitted
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
			% fixed population value. The length is therefore capped at 20
			% samples (and at least 10), so the cap only affects longer series.
			l = max(min(20, floor(N * samplep(2))), 10);
		else % specified an absolute interval
			l = samplep(2);
		end

		% numPred random starting points (uniform on 1..N-l+1), reproducible from the seed:
		spts = 1 + floor((N - l + 1) * BF_Random(numPred, BF_RandomSeed(randomSeed)));
		r(:, 1) = spts;
		r(:, 2) = spts + l - 1;

	case 'spread'
		if samplep(2) < 1 % specified a fraction of time series
			l = max(min(20, floor(N * samplep(2))), 10); % (capped, as for 'rand')
		else % specified an absolute interval
			l = samplep(2);
		end
		% start points on 1..N-l+1 from the golden-ratio (Weyl) sequence frac(j*phi):
		spts = 1 + floor((N - l + 1) * mod((1:numPred)' * 0.6180339887498949, 1));
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

	% Step-ahead predictions across the test segment, using the model fitted to all the
	% data and a diffuse initial state (the first predictions are undefined):
	yp = BF_SSPredict(A, K, C, yTest, steps);
	defined = ~isnan(yp);
	if sum(defined) < 3
		rmserrs(i) = NaN; ac1s(i) = NaN; meandiffs(i) = NaN; stdrats(i) = NaN;
		continue
	end
	yp = yp(defined); yTest = yTest(defined);

	% Get statistics on residuals:
	mres = yp - yTest;

	rmserrs(i) = sqrt(mean(mres.^2));
	if ~(std(mres) > 1e-9 * max(sigmaY, mean(abs(mres)))) % (the rule of MF_ResidualAnalysis)
		% errors that do not vary beyond the level of rounding (e.g., the segment is
		% predicted exactly): their autocorrelation is not meaningful, and the segment is
		% left out of the ac1 summaries
		ac1s(i) = NaN;
	else
		ac1s(i) = CO_AutoCorr(mres, 1, 'Fourier');
	end

	% Get statistics on output time series
	meandiffs(i) = abs(mean(yp) - mean(yTest));
	% Guard against near-constant test segments: std(yTest) can land at
	% floating-point noise (rather than a genuinely small but real value)
	% when a short random/uniform segment happens to sit in a flat/clipped
	% stretch of the input. Then std(yp)/std(yTest) is undefined rather
	% than just large, so exclude it instead of letting it dominate the
	% mean/std/iqr summaries below.
	if std(yTest) < 1e-6 * sigmaY
		stdrats(i) = NaN;
	else
		stdrats(i) = std(yp) / std(yTest);
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
% (The median is not returned, being largely redundant with the mean; the mean
%  absolute error, mabserrs, is not returned at all, since it adds no dimension
%  beyond the RMS error.)
out.stde_mean = mean(rmserrs);
out.stde_std = std(rmserrs);
out.stde_iqr = iqr(rmserrs);

% Autocorrelations at lag 1, ac1s
% NOT absolute values of ac1s... absolute values of operations on *raw* ac1s...
validAC1s = ac1s(~isnan(ac1s));
out.ac1_mean = abs(mean(validAC1s));
out.ac1_median = abs(median(validAC1s));
out.ac1_std = std(validAC1s);
out.ac1_iqr = iqr(validAC1s);

% Differences in mean between two series
% (The median is not returned, being largely redundant with the mean.)
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
