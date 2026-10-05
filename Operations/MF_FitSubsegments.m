function out = MF_FitSubsegments(y, model, order, subsetHow, samplep, randomSeed)
% MF_FitSubsegments   How the fit of a model varies across segments of the time series.
%
% Fits the same kind of model to many segments of the time series, and summarizes how
% the fitted parameters and goodness of fit vary from segment to segment. The spread
% of the parameters (including in-sample goodness-of-fit statistics) indicates
% stationarity, and the values of goodness of fit indicate the suitability of the
% model. This code inherits strongly from MF_CompareTestSets.
%
% ---INPUTS:
% y, the input time series.
%
% model, the model to fit in each segment of the time series:
%           'arsbc': fits an AR model of the best order (1 to 10) by the Schwarz
%                       Bayesian criterion (SBC) using the ARfit package. Outputs
%                       are on how the optimal order and the SBC vary in different
%                       parts of the time series. (The order input is not used.)
%           'ar': fits an AR model of a specified order using ar from MATLAB's
%                   System Identification Toolbox. Outputs are on how Akaike's Final
%                   Prediction Error (FPE) and the fitted AR parameters vary across
%                   the segments.
%           'ss': fits a state-space model of a given order using n4sid from
%                   MATLAB's System Identification Toolbox. Outputs are on how the
%                   FPE varies.
%           'arma': fits an ARMA model using armax from MATLAB's System
%                   Identification Toolbox. Outputs are on the FPE and the fitted AR
%                   and MA parameters.
%           'arcrosspred': splits the series into samplep non-overlapping segments,
%                   fits an AR model of the given order to each, and uses every
%                   segment's model to predict, one step ahead, every segment
%                   (including itself). This forms a samplep x samplep matrix of
%                   cross-prediction root-mean-square errors, a linear-model
%                   analogue of SY_nstat_z's nonlinear zeroth-order cross-prediction
%                   matrix. Only the spread and off-diagonal statistics of that
%                   matrix are returned. Requires subsetHow = 'uniform' and a scalar
%                   samplep (a segment count).
%           The default is 'ss'.
%
% order, the order of the model to fit (default 2; a two-element vector [p, q] for
%           'arma').
%
% subsetHow, how to choose segments from the time series, either 'uniform'
%           (evenly spaced) or 'rand' (at random) (default).
%
% samplep, a two-vector specifying how many segments to take and of what length, of
%           the form [nsamples, length], where length can be a proportion of the
%           time-series length (default [20, 0.1], i.e., 20 segments of 10% of the
%           time-series length). For model = 'arcrosspred', a scalar giving the
%           number of non-overlapping segments to partition the whole series into.
%
% randomSeed, the seed of the random start points (see BF_RandomSeed; the numbers come
%           from the portable generator BF_Random; for when subsetHow is 'rand')
%
% ---OUTPUTS: depend on the model.
% For 'arsbc', statistics across segments of the best AR order and its SBC:
%   orders_mode, orders_mean, orders_std, orders_max, orders_min, orders_range:
%       the mode, mean, standard deviation, maximum, minimum and range of the order
%   sbcs_mean, sbcs_std, sbcs_range, sbcs_min, sbcs_max: the mean, standard
%       deviation, range, minimum and maximum of the SBC of the best-order fit
% For 'ar', 'ss' and 'arma', statistics across segments of the FPE:
%   fpe_std, fpe_mean, fpe_max, fpe_min, fpe_range: the standard deviation, mean,
%       maximum, minimum and range of the FPE
% For 'ar', statistics across segments of the fitted AR coefficients, as in the
% polynomial 1 + a_1 z^-1 + ... (the negative of the usual AR coefficients), for each
% lag k up to the order:
%   a_k_std, a_k_mean, a_k_max, a_k_min (e.g., a_1_std, a_1_mean, a_1_max, a_1_min,
%       a_2_std, a_2_mean, a_2_max, a_2_min for order 2)
% For 'arma', the same statistics of the AR coefficients, p_k_std, p_k_mean, p_k_max,
% p_k_min (k = 1 to order(1)), and of the MA coefficients, q_k_std, q_k_mean, q_k_max,
% q_k_min (k = 1 to order(2)).
% For 'arcrosspred', statistics of the cross-prediction error matrix (a row for each
% predicting model, a column for each predicted segment):
%   std, range, iqr: the standard deviation, range and interquartile range over all
%       entries
%   stdoffdiag, rangeoffdiag, iqroffdiag: the same over the off-diagonal entries
%   stdmean, rangemean, stdmedian, rangemedian: the standard deviation and range,
%       across predicted segments, of the mean and of the median error
%   rangerange, stdrange, rangestd, stdstd: the range or standard deviation, across
%       predicted segments, of the range or of the standard deviation of the errors
%       (rangerange: range of range, stdrange: std of range, rangestd: range of std,
%       stdstd: std of std)
%   mineig: the smallest real part of the eigenvalues of the matrix
%
% ---NOTES:
% 'arcrosspred': the level statistics (trace, mean, min, max, eigenvalue levels) of
% the cross-prediction matrix correlated at |r| = 0.89-0.98 with this function's own
% 'ar' fpe_mean/min/max fields on two collections of real-world series (150 series
% each), so are not returned. The spread and off-diagonal statistics, which correlated
% only at |r| = 0.4-0.88 with fpe_* and a_1_*, and 0.4-0.85 with SY_nstat_z's spread
% statistics, are kept.
%
% The 'arma' registration (order = [2, 2], 25 uniform 10%-length segments) was
% deregistered: its AR-driven fields correlated at |r| = 0.69-0.99 with
% the much cheaper 'ar' registration, and the MA-driven fields (q_k_*) were noisy.
% The 'ar', 'arsbc' and 'ss' registrations are unaffected.

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
%% Check Inputs
% ------------------------------------------------------------------------------

% (1) y: column vector time series
if nargin < 1 || isempty(y)
	error('Give us a time series, ya mug');
end
% Convert y to time series object
y = iddata(y, [], 1);

% (2) model, the type of model to fit
if nargin < 2 || isempty(model)
	model = 'ss'; % fit a state space model by default
end

% (3) order of model, order
if nargin < 3 || isempty(order)
	order = 2; % model of order 2 by default. Not very good defaults.
end

% (4) How to choose subsets from the time series, subsetHow
if nargin < 4 || isempty(subsetHow)
	subsetHow = 'rand'; % takes segments randomly from time series
end

% (5) Sampling parameters, samplep
if nargin < 5 || isempty(samplep)
	samplep = [20, 0.1]; % sample 20 times with 10%-length subsegments
end

% (6) randomSeed: how to treat the randomization
if nargin < 6
	randomSeed = [];
end

% 'arcrosspred' needs a genuine non-overlapping partition of the series
% (not resampled/overlapping segments) for cross-prediction to make sense:
if strcmp(model, 'arcrosspred') && (~strcmp(subsetHow, 'uniform') || numel(samplep) ~= 1)
	error('''arcrosspred'' requires subsetHow = ''uniform'' and a scalar samplep (a non-overlapping segment count)');
end

% ------------------------------------------------------------------------------
%% Set the ranges beforehand
% ------------------------------------------------------------------------------
% Number of samples to take, numPred
numPred = samplep(1);
r = zeros(numPred, 2); % ranges

switch subsetHow
	case 'rand'
		if samplep(2) < 1 % specified a fraction of time series
			l = floor(N * samplep(2));
		else % specified an absolute interval
			l = samplep(2);
		end

		% numPred random starting points (uniform on 1..N-l+1), reproducible from the seed:
		spts = 1 + floor((N - l + 1) * BF_Random(numPred, BF_RandomSeed(randomSeed)));
		r(:, 1) = spts;
		r(:, 2) = spts + l - 1;

	case 'uniform'
		if length(samplep) == 1 % size will depend on number of unique subsegments
			spts = round(linspace(0, N, numPred + 1)); % numPred+1 boundaries = numPred portions
			r(:, 1) = spts(1:numPred) + 1;
			r(:, 2) = spts(2:end);
		else
			if samplep(2) < 1 % specified a fraction of time series
				l = floor(N * samplep(2));
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

% ------------------------------------------------------------------------------
%% Fit the model to each training set
% ------------------------------------------------------------------------------
% model will be stored as a model object, m
% model is fitted using the entire dataset as the training set
% test sets will be smaller chunks of this.
% [could also fit multiple models using data not in multiple test sets, but
% this is messier]
switch model
	case 'arsbc'
		%% Fit AR models of 'best' order according to SBC, using arfit package
		% fit AR models of 'best' order, return statistics on how this best
		% order changes. The order input argument is not used for this
		% option.
		orders = zeros(numPred, 1);
		sbcs = zeros(numPred, 1);
		yy = y.y;
		for i = 1:numPred
			% Use arfit software to retrieve the optimum AR(p) order by
			% Schwartz's Bayesian Criterion, SBC (or BIC); in the range
			% p = 1-10
			% Enforce zero mean level. This could be relaxed.
			try
				[west, Aest, Cest, SBC] = ARFIT_arfit(yy(r(i, 1):r(i, 2)), 1, 10, 'sbc', 'zero');
			catch emsg
				if contains(emsg.message, 'too short') % (ARFIT_arfit: 'Time series (N = %u) too short.')
					fprintf(1, 'Time Series is too short for ARFIT\n');
					out = NaN; return
				else
					error('Problem fitting AR model');
				end
			end
			orders(i) = length(Aest);
			sbcs(i) = min(SBC);
		end

		% Return statistics
		out.orders_mode = mode(orders);
		out.orders_mean = mean(orders);
		out.orders_std = std(orders);
		out.orders_max = max(orders);
		out.orders_min = min(orders);
		out.orders_range = range(orders);

		out.sbcs_mean = mean(sbcs);
		out.sbcs_std = std(sbcs);
		out.sbcs_range = range(sbcs);
		out.sbcs_min = min(sbcs);
		out.sbcs_max = max(sbcs);

	case 'ar'
		%% Fit AR model of specified order
		% Return statistics on parameters and goodness of fit

		%% Check that a System Identification Toolbox license is available to run the 'ar' function:
		BF_CheckToolbox('identification_toolbox');

		fpes = zeros(numPred, 1);
		as = zeros(numPred, order + 1);
		for i = 1:numPred
			% fit the ar model
			m = ar(y(r(i, 1):r(i, 2)), order);
			% get parameters and goodness of fit
			fpes(i) = m.EstimationInfo.FPE;
			as(i, :) = m.a;
		end

		% statistics on FPE
		out.fpe_std = std(fpes);
		out.fpe_mean = mean(fpes);
		out.fpe_max = max(fpes);
		out.fpe_min = min(fpes);
		out.fpe_range = range(fpes);

		% Statistics on fitted AR parameters
		for i = 1:order % first column will be ones
			% Dynamic field referencing:
			out.(['a_', num2str(i), '_std']) = std(as(:, i + 1));
			out.(['a_', num2str(i), '_mean']) = mean(as(:, i + 1));
			out.(['a_', num2str(i), '_max']) = max(as(:, i + 1));
			out.(['a_', num2str(i), '_min']) = min(as(:, i + 1));
			% eval(sprintf('out.a_%u_std = std(as(:,%u+1));',i,i));
			% eval(sprintf('out.a_%u_mean = mean(as(:,%u+1));',i,i));
			% eval(sprintf('out.a_%u_max = max(as(:,%u+1));',i,i));
			% eval(sprintf('out.a_%u_min = min(as(:,%u+1));',i,i));
		end

	case 'arcrosspred'
		%% Fit an AR model to each of numPred non-overlapping segments, then
		%% use every segment's model to 1-step-ahead predict every segment
		%% (including itself), forming a numPred x numPred cross-prediction
		%% RMSE matrix. Only its spread/off-diagonal statistics are
		%% returned -- see NOTES for why the level statistics are omitted.

		%% Check that a System Identification Toolbox license is available:
		BF_CheckToolbox('identification_toolbox');

		minSegLength = 5 * (order + 1); % need enough points to fit AR(order) and predict meaningfully
		if any(r(:, 2) - r(:, 1) + 1 < minSegLength)
			warning('Segments too short to reliably cross-predict with an AR(%u) model', order);
			out = NaN; return
		end

		try
			segs = cell(numPred, 1);
			models = cell(numPred, 1);
			for i = 1:numPred
				segs{i} = y(r(i, 1):r(i, 2));
				models{i} = ar(segs{i}, order);
			end
			xperr = zeros(numPred); % cross-prediction RMSE from using segment i's model on segment j
			for i = 1:numPred
				for j = 1:numPred
					yp = predict(models{i}, segs{j}, 1);
					res = yp.y - segs{j}.y;
					xperr(i, j) = sqrt(mean(res.^2));
				end
			end
		catch
			% A segment was degenerate (e.g. near-constant) for AR fitting/prediction
			out = NaN; return
		end

		out.std = std(xperr(:));
		out.range = range(xperr(:));
		out.iqr = iqr(xperr(:));

		lowertri = tril(xperr, -1); lowertri = lowertri(lowertri > 0);
		uppertri = triu(xperr, 1); uppertri = uppertri(uppertri > 0);
		offdiag = [lowertri; uppertri];
		if isempty(offdiag)
			out.iqroffdiag = NaN;
			out.stdoffdiag = NaN;
			out.rangeoffdiag = NaN;
		else
			out.iqroffdiag = iqr(offdiag);
			out.stdoffdiag = std(offdiag);
			out.rangeoffdiag = range(offdiag);
		end

		% Comparing columns/rows (i.e., how differently each segment's model
		% behaves as a predictor vs. as a target)
		out.stdmean = std(mean(xperr));
		out.rangemean = range(mean(xperr));
		out.stdmedian = std(median(xperr));
		out.rangemedian = range(median(xperr));
		out.rangerange = range(range(xperr));
		out.stdrange = std(range(xperr));
		out.rangestd = range(std(xperr));
		out.stdstd = std(std(xperr));

		% Eigenvalues
		realEigs = real(eig(xperr));
		out.mineig = min(realEigs);

	case 'ss'
		%% Fit state space models of specified order
		% Return statistics on goodness of fit
		% Could do parameters too, but I this would involve many outputs
		fpes = zeros(numPred, 1);
		for i = 1:numPred
			try  m = n4sid(y(r(i, 1):r(i, 2)), order);
			catch
				% Some range of the time series is invalid for fitting the
				% model to.
				error('Couldn''t fit this state space model')
			end
			fpes(i) = m.EstimationInfo.FPE;
		end

		% statistics on FPE
		out.fpe_std = std(fpes);
		out.fpe_mean = mean(fpes);
		out.fpe_max = max(fpes);
		out.fpe_min = min(fpes);
		out.fpe_range = range(fpes);

	case 'arma'
		%% fit an ARMA model of specified order(s)
		% Note: order should be a two-component vector
		% Output parameters and goodness of fit
		fpes = zeros(numPred, 1);
		ps = zeros(numPred, order(1) + 1);
		qs = zeros(numPred, order(2) + 1);

		for i = 1:numPred
			try
				m = armax(y(r(i, 1):r(i, 2)), order);
			catch emsg
				error('Couldn''t fit this ARMA model')
			end
			fpes(i) = m.EstimationInfo.FPE;
			ps(i, :) = m.a;
			qs(i, :) = m.c;
		end

		% statistics on FPE
		out.fpe_std = std(fpes);
		out.fpe_mean = mean(fpes);
		out.fpe_max = max(fpes);
		out.fpe_min = min(fpes);
		out.fpe_range = range(fpes);

		% Statistics on fitted AR parameters, p
		for i = 1:order(1) % first column will be ones
			out.(['p_', num2str(i), '_std']) = std(ps(:, i + 1));
			out.(['p_', num2str(i), '_mean']) = mean(ps(:, i + 1));
			out.(['p_', num2str(i), '_max']) = max(ps(:, i + 1));
			out.(['p_', num2str(i), '_min']) = min(ps(:, i + 1));
		end

		% Statistics on fitted MA parameters, q
		for i = 1:order(2) % first column will be ones
			out.(['q_', num2str(i), '_std']) = std(qs(:, i + 1));
			out.(['q_', num2str(i), '_mean']) = mean(qs(:, i + 1));
			out.(['q_', num2str(i), '_max']) = max(qs(:, i + 1));
			out.(['q_', num2str(i), '_min']) = min(qs(:, i + 1));
		end
	otherwise
		error('Unknown model ''%s''', model);
end

end
