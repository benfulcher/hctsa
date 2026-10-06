function out = MF_GP_LocalPrediction(y, covFunc, numTrain, numTest, numPreds, pmode, randomSeed, numSplits)
% MF_GP_LocalPrediction   How well a Gaussian process fitted to short windows of the series predicts nearby held-out values.
%
% Takes numPreds windows spread evenly along the time series. In each window, a
% Gaussian process (GP) with the covariance function covFunc is fitted to a
% training segment (with the noise standard deviation bounded below by 1% of that of
% the training data) and used to predict numTest held-out values, which are compared
% with the true values. Each window is first standardized using the mean and
% standard deviation of its training data. The outputs summarize the prediction
% errors (absolute, and relative to the GP's own 95% error bar), the size of the
% error bars, the fitted hyperparameters across windows, and the per-point
% negative log marginal likelihoods.
%
% Uses GP fitting code from the gpml toolbox, which is available here:
% http://gaussianprocess.org/gpml/code.
%
% ---INPUTS:
% y, the input time series
%
% covFunc, covariance function in the standard form for the gpml package.
%           E.g., covFunc = {'covSum', {'covSEiso','covNoise'}} combines squared
%           exponential and noise terms (the default)
%
% numTrain, the number of training samples (for each window; default: 20)
%
% numTest, the number of test samples (for each window; default: 5)
%
% numPreds, the number of windows, and so of predictions made (default: 10)
%
% pmode, the prediction mode (default: 'frombefore'):
%       (i) 'beforeafter': trains on numTrain samples on each side of a gap of
%                           numTest samples, and predicts the samples in the gap,
%       (ii) 'frombefore': trains on numTrain samples and predicts the numTest
%                    samples that follow,
%       (iii) 'randomgap': trains on a random numTrain of the numTrain + numTest
%                    samples in the window and predicts the other numTest samples, and
%       (iv) 'spreadgap': as 'randomgap', but with deterministic splits: the test
%                    sets are taken, in the evenly spread order of BF_SpreadPerm, from
%                    the list of all nchoosek(numTrain + numTest, numTest) possible test
%                    sets (in lexicographic order), a different one for each fit (cycling
%                    through the list when there are more fits than test sets). The
%                    fits cover the possible splits evenly, as random splits do on
%                    average, and the outputs do not depend on a seed.
%
% randomSeed, the seed of the random splits (see BF_RandomSeed; they come from the
%               portable generator BF_Random; for 'randomgap' prediction)
%
% numSplits, the number of different splits made of each window, for 'randomgap' and
%               'spreadgap' prediction (default: 8). Each split is fitted and predicted as
%               a window of its own, so the outputs summarize numPreds * numSplits fits: a
%               few splits of ten windows are dominated by which splits were made.
%
% ---OUTPUTS:
% meanabs_run, maxabs_run, minabs_run: mean, maximum, and minimum over windows of
%       the mean absolute prediction error in a window (in units of the training
%       data's standard deviation)
% meanabs_std_run, maxabs_std_run, minabs_std_run: the same, with each error in
%       units of the GP's 95% error bar (twice its predictive standard deviation)
% q90abs_run, q90abs_std_run: the 90th percentile (MATLAB's quantile) over windows of
%       the mean absolute prediction error in a window, without and with each error
%       in units of the 95% error bar
% q10abs_run: the 10th percentile over windows of the mean absolute prediction error
%       in a window
% low25abs_std_run: the mean over the lowest quarter of the windows (the ceil(n/4)
%       smallest of the n values) of the mean absolute prediction error in a window,
%       in units of the 95% error bar
%       (These four summarize the upper and lower tails robustly: the maximum and
%       minimum over windows are each set by a single fit, often a badly conditioned
%       one, and so mostly reflect which splits were made.)
% meanabs, maxabs, minabs: mean, maximum, and minimum over all predicted points of
%       the absolute prediction error
% meanabs_std, maxabs_std, minabs_std: the same, in units of the 95% error bar
% maxerrbar, meanerrbar, minerrbar: maximum, mean, and minimum over all predicted
%       points of the 95% error bar half-width (twice the predictive standard
%       deviation)
% high25errbar: the mean of the largest quarter of these error bars (over all
%       predicted points)
% meanlogh1, meanlogh2, meanlogh3, stdlogh1, stdlogh2, stdlogh3: mean and standard
%       deviation across windows of each log hyperparameter of the fitted
%       covariance function (for squared exponential plus noise: log length scale,
%       log amplitude, log noise standard deviation)
% maxnlml, minnlml, stdnlml: maximum, minimum, and standard deviation across
%       windows of the negative log marginal likelihood of the fitted model on the
%       training data of the window, divided by the number of training points
% q90nlml: the 90th percentile of the same across windows
%
% For 'randomgap' and 'spreadgap', each split counts as a window in these summaries.
%
% ---NOTES:
% Windows whose training data are constant (standard deviation below 1e-8 of that of
% the series) cannot be standardized, and are left out of every statistic; the output
% is NaN if every window is left out.
% The 'standard errors' in the code (stderrs) are 2*sqrt(S2), i.e., 95% error
% bars, so the outputs ending in _std are in units of these, not of one standard
% deviation. The predictive variance S2 includes the likelihood noise.
% For 'randomgap', all the random splits are drawn once, from one seed, before the
% loop over windows: the sequence of splits is the same on every run (with a fixed seed).
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
doPlot = false; % whether to plot outputs
N = length(y); % time-series length

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if size(y, 2) > size(y, 1)
	y = y'; % ensure a column vector input
end
if nargin < 2 || isempty(covFunc),
	fprintf(1, 'Using a default covariance function: sum of squared exponential and noise\n');
	covFunc = {'covSum', {'covSEiso', 'covNoise'}};
end

if nargin < 3 || isempty(numTrain)
	numTrain = 20; % 20 previous data points to predict the next
end

if nargin < 4 || isempty(numTest)
	numTest = 5; % test on 5 data points into the future
end

if nargin < 5 || isempty(numPreds) % number of predictions
	numPreds = 10; % do it 10 times (equally-spaced) through the time series
end

if nargin < 6 || isempty(pmode)
	pmode = 'frombefore'; % predicts from previous numTrain datapoints
	% can also be 'randomgap' -- fills in random gaps in the middle of a string of
	% data
end

% randomSeed: how to treat the randomization
if nargin < 7
	randomSeed = [];
end

% numSplits: splits of each window ('randomgap' and 'spreadgap' only)
if nargin < 8 || isempty(numSplits)
	numSplits = 8;
end
if ~ismember(pmode, {'randomgap', 'spreadgap'})
	numSplits = 1;
end

% ------------------------------------------------------------------------------
%% Set up loop
% ------------------------------------------------------------------------------
if ismember(pmode, {'frombefore', 'randomgap', 'spreadgap'})
	spns = floor(linspace(1, N - (numTest + numTrain), numPreds)); % starting positions
elseif strcmp(pmode, 'beforeafter')
	spns = floor(linspace(1, N - (numTest + numTrain * 2), numPreds)); % starting positions
end
% Each window is used numSplits times (with a different split each time):
spns = repelem(spns, numSplits);
numPreds = numPreds * numSplits;

% Details of GP:
meanFunc = {'meanZero'}; % zero-mean process
likFunc = @likGauss; % likelihood function (Gaussian)
% Exact Gaussian inference. (Was @infLaplace: with a Gaussian likelihood the
% Laplace approximation is exact but is computed by Newton iteration, and
% gpml's infLaplace warm-starts that iteration from a PERSISTENT copy of the
% previous call's solution, so every fit depended on whatever series was
% fitted before it -- the same series gave outputs differing in the 4th digit
% from call to call once the optimizer amplified the difference.)
infAlg = @infGaussLik;
nfevals = -50;
hyp = struct; % structure for storing hyperparameter information in latest version of GMPL toolbox

% Initialize variables:
mus = zeros(numTest, numPreds); % predicted values
stderrs = zeros(numTest, numPreds); % standard errors on predictions
yss = zeros(numTest, numPreds); % test values
nlmls = NaN(numPreds, 1); % negative log marginal likelihoods of model, per training point

nhps = eval(feval(covFunc{:})); % number of hyperparameters
loghypers = zeros(nhps, numPreds); % loghyperparameters

% The train/test splits, one column per fit (so that the windows get different splits):
if strcmp(pmode, 'spreadgap')
	% Deterministic splits: the test sets are taken from the list of all
	% nchoosek(numTrain + numTest, numTest) possible ones (in lexicographic order) in the
	% evenly spread order of BF_SpreadPerm, a different test set for each of the
	% numPreds (windows x splits) fits (cycling through the list if it is shorter):
	numSplitsAll = nchoosek(numTrain + numTest, numTest);
	if numSplitsAll > 1e6
		error('Too many possible splits (%g) for ''spreadgap''', numSplitsAll);
	end
	testSets = nchoosek(1:numTrain + numTest, numTest);
	spreadOrder = BF_SpreadPerm(numSplitsAll);
	splits = zeros(numTrain + numTest, numPreds);
	for i = 1:numPreds
		testSet = testSets(spreadOrder(mod(i - 1, numSplitsAll) + 1), :);
		splits(:, i) = [setdiff(1:numTrain + numTest, testSet), testSet]';
	end
elseif strcmp(pmode, 'randomgap')
	% Random splits, reproducible from the seed:
	[~, splits] = sort(reshape(BF_Random((numTrain + numTest) * numPreds, BF_RandomSeed(randomSeed)), ...
										numTrain + numTest, numPreds));
end

for i = 1:numPreds
	%% (0) Set up test and training sets
	switch pmode
		case 'frombefore'
			tt = (1:numTrain)'; % times (make from 1)
			rt = spns(i):spns(i) + numTrain - 1; % training range
			yt = y(rt); % training data

			ts = (numTrain + 1:numTrain + 1 + numTest - 1)'; % times
			rs = spns(i) + numTrain:spns(i) + numTrain + numTest - 1; % test range
			ys = y(rs); % test data

		case {'randomgap', 'spreadgap'}
			t = (1:numTrain + numTest)';
			r = splits(:, i)';
			yy = y(spns(i):spns(i) + numTrain + numTest - 1);

			rt = sort(r(1:numTrain), 'ascend');
			tt = t(rt);
			yt = yy(rt);

			rs = sort(r(numTrain + 1:end), 'ascend');
			ts = t(rs);
			ys = yy(rs);

		case 'beforeafter'
			t = (1:2 * numTrain + numTest)';
			yy = y(spns(i):spns(i) + 2 * numTrain + numTest - 1);

			rt = [1:numTrain, numTrain + numTest + 1:numTrain * 2 + numTest];
			tt = t(rt);
			yt = yy(rt);

			rs = (numTrain + 1:numTrain + numTest);
			ts = t(rs);
			ys = yy(rs);

		otherwise
			error('Unknown prediction mode ''%s''', pmode);
	end

	% A window whose training data are constant (to rounding error, relative to the
	% series) cannot be standardized and carries no information about a GP: skip it
	if ~(std(yt) > 1e-8 * std(y))
		continue
	end

	% Process to normalize scales
	ys = (ys - mean(yt)) / std(yt); % same transformation as training set
	yt = (yt - mean(yt)) / std(yt); % zscore training set

	% ------------------------------------------------------------------------------
	%% (1) Learn hyperparameters from training set (t)
	% ------------------------------------------------------------------------------

	% Initialize mean and likelihood
	hyp.mean = []; hyp.lik = log(0.1);
	hyp.cov = [];

	% loghyper = MF_GP_LearnHyperp(covFunc,-50,tt,yt);
	try
		hyp = MF_GP_LearnHyperp(tt, yt, covFunc, meanFunc, likFunc, infAlg, nfevals, hyp);
	catch emsg
		fprintf(1, 'Unable to learn hyperparameters for this time series\n');
		out = NaN; return
	end
	if ~isstruct(hyp) % MF_GP_LearnHyperp returns NaN (not a struct) when the data isn't suited to GP fitting
		fprintf(1, 'Unable to learn hyperparameters for this time series\n');
		out = NaN; return
	end
	loghyper = hyp.cov;

	if any(isnan(loghyper))
		fprintf(1, 'Unable to learn hyperparameters for this time series\n');
		out = NaN; return
	end

	loghypers(:, i) = loghyper;

	% Get marginal likelihood for this model with hyperparameters optimized
	% over training data
	% nlmls(i) = gpr(loghyper, covFunc, tt, yt);
	% (negative log marginal likelihood, gpml's nlZ, divided by the number of training points)
	nlmls(i) = gp(hyp, infAlg, meanFunc, covFunc, likFunc, tt, yt) / length(tt);

	% ------------------------------------------------------------------------------
	%% (2) Evaluate at test set (s)
	% ------------------------------------------------------------------------------

	% Evaluate at test points based on training time/data, predicting for
	% test times/data
	% [mu, S2] = gpr(loghyper, covFunc, tt, yt, ts); % old version
	[mu, S2] = gp(hyp, infAlg, meanFunc, covFunc, likFunc, tt, yt, ts); % evaluate at new time points, ts

	% Compare to actual test data --> store in row of errs
	mus(:, i) = mu; % ~predicted values for time series points
	stderrs(:, i) = 2 * sqrt(S2); % ~errors on those predictions
	yss(:, i) = ys;

	% Plot
	if doPlot
		if strcmp(pmode, 'frombefore')
			plot(tt, yt, '.-k');
			hold on;
			plot(ts, ys, '.-b');
			errorbar(ts, mu, 2 * sqrt(S2), 'm');
			hold off;
		else
			plot(tt, yt, 'ok');
			hold on;
			plot(ts, ys, 'ob');
			errorbar(ts, mu, 2 * sqrt(S2), 'm');
			hold off;
		end
	end

	%     for j=1:numTest
	%         % set up structure output
	%         err = abs(mu(j)-ys(j))/sqrt(S2(j)); % in units of std at this point
	%         eval(['out.abserr' num2str(i) '_' num2str(j) ' = err;']);
	%     end

end

% Ok, we're done.

% Drop the skipped windows (those with constant training data)
keep = ~isnan(nlmls);
if ~any(keep)
	out = NaN; return
end
mus = mus(:, keep); stderrs = stderrs(:, keep); yss = yss(:, keep);
loghypers = loghypers(:, keep); nlmls = nlmls(keep);

% ------------------------------------------------------------------------------
%% Return statistics on how well it did
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% (1) PREDICTION ERROR MEASURES
% ------------------------------------------------------------------------------

% Absolute errors:
allabserrs = abs(mus - yss);
% In units of standard errors (95% confidence interval error bars)
allstderrs = abs(mus - yss) ./ stderrs;

% ---
% * Stats on all errors:
% ---

% largest error:
out.maxabs_std = max(allstderrs(:));
out.maxabs = max(allabserrs(:));

% smallest error:
out.minabs_std = min(allstderrs(:));
out.minabs = min(allabserrs(:));

% mean error (across all):
out.meanabs_std = mean(allstderrs(:));
out.meanabs = mean(allabserrs(:));

% ---
% * Stats on errors per run
% ---

% Summary of how it did on each run:
stderr_run = mean(allstderrs);
abserr_run = mean(allabserrs);

% Mean error for a run
out.meanabs_std_run = mean(stderr_run);
out.meanabs_run = mean(abserr_run);

% Max error for a run
out.maxabs_std_run = max(stderr_run);
out.maxabs_run = max(abserr_run);

% Min error for a run
out.minabs_std_run = min(stderr_run);
out.minabs_run = min(abserr_run);

% Upper (and lower) tails over the fits: the maximum (minimum) over the fits is set by
% one fit (often a single badly conditioned one), and so depends on which splits were
% made; a high (low) quantile, or the mean of the highest (lowest) quarter, keeps the
% meaning ('how large/small does it get') and is reproducible
out.q90abs_run = quantile(abserr_run, 0.9);
out.q10abs_run = quantile(abserr_run, 0.1);
out.q90abs_std_run = quantile(stderr_run, 0.9);
out.low25abs_std_run = SUB_tailMean(stderr_run, 'low');

% Error bar stats:
out.maxerrbar = max(stderrs(:)); % largest error bar
out.high25errbar = SUB_tailMean(stderrs(:), 'high'); % mean of the largest quarter of the error bars
out.meanerrbar = mean(stderrs(:)); % mean error bar length
out.minerrbar = min(stderrs(:)); % minimum error bar length

% ------------------------------------------------------------------------------
%% (2) HYPERPARAMETER MEASURES
% ------------------------------------------------------------------------------
% mean and std for each hyperparameter
for i = 1:nhps
	out.(sprintf('meanlogh%u', i)) = mean(loghypers(i, :));
	out.(sprintf('stdlogh%u', i)) = std(loghypers(i, :));
end

% ------------------------------------------------------------------------------
%% (3) Marginal likelihood measures
% ------------------------------------------------------------------------------
% Worst (maximum) per-point marginal neg-log-likelihood attained
% Best (minimum) per-point marginal neg-log-likelihood attained
% spread in per-point marginal neg-log-likelihoods
% (Previously maxmlik, minmlik and stdmlik: the log marginal likelihood, i.e., the
%  negative of nlZ, summed over the training points rather than per point.)

out.maxnlml = max(nlmls);
out.q90nlml = quantile(nlmls, 0.9);
out.minnlml = min(nlmls);
out.stdnlml = std(nlmls);

% ------------------------------------------------------------------------------
function m = SUB_tailMean(x, whichTail)
	% Mean of the highest ('high') or lowest ('low') quarter of the values in x
	% (the ceil(n/4) largest or smallest of the n values)
	x = sort(x(:), 'ascend');
	k = ceil(length(x) / 4);
	if strcmp(whichTail, 'high')
		m = mean(x(end - k + 1:end));
	else
		m = mean(x(1:k));
	end
end
% ------------------------------------------------------------------------------

end
