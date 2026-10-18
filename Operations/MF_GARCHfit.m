function out = MF_GARCHfit(y, preproc, P, Q, randomSeed, modelType, innovationDist)
% MF_GARCHfit   A GARCH model of the changing variance of the series, and what it leaves in the residuals.
%
% Fits a Generalized Autoregressive Conditional Heteroskedasticity (GARCH) model to a
% time series, namely:
%
% (1) Preprocessing the data to remove strong trends (and, optionally, to whiten it),
% (2) Pre-estimation to calculate initial correlation properties of the time
%       series and motivate a GARCH model,
% (3) Fitting a GARCH model, returning goodness-of-fit statistics and parameters
%           of the fitted model, and
% (4) Post-estimation, involving statistics on the residuals and the standardized
%           residuals.
%
% The preprocessing is BF_Whiten: for preproc = 'ar', the series is detrended and
% replaced by the preprocessing (from PP_PreProcess) that maximizes whiteness under an
% AR(2) model, if one improves on the detrended series by more than 5%; for preproc =
% 'none' the series is left as it is (not detrended). The result is then z-scored, so
% the model is of the variance around a zero mean.
%
% The model is fitted by maximum likelihood with BF_GARCHFit, which specifies the
% likelihood (Bollerslev's GARCH, the threshold form of Glosten, Jagannathan and
% Runkle, and Student-t innovations), its bounds, the presample values and a
% deterministic, converged optimizer; no toolbox is needed for the fit.
%
% ---INPUTS:
% y, the input time series
%
% preproc, the preprocessing to apply, 'ar' (default) or 'none'
%
% P, the GARCH model order (the number of lagged variances; default 1)
%
% Q, the ARCH model order (the number of lagged squared innovations; default 1)
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (for pre-processing: PP_PreProcess)
%
% modelType, the conditional variance model to fit: 'garch' (default,
%               symmetric response to shocks) or 'gjr' (GJR-GARCH, adds a
%               leverage/asymmetry term so negative and positive shocks can
%               have different effects on variance).
%
% innovationDist, the assumed innovation distribution: 'gaussian' (default)
%               or 't' (Student's t, estimates a degrees-of-freedom
%               parameter to capture fat tails beyond what GARCH-filtering
%               alone accounts for).
%
% ---OUTPUTS:
% hasStructure: 1 if the fitted variance model beats constant variance in the
%       likelihood-ratio test described below, otherwise 0
% Fitted parameters (standard errors have 'err' in the name):
% constant, constanterr: the constant term of the variance equation, and its error
% offset: the mean offset of the model (0 for a z-scored series)
% GARCH_1, GARCHerr_1, ..., GARCH_P, GARCHerr_P: the coefficients of the P lagged
%       variances, and their errors
% ARCH_1, ARCHerr_1, ..., ARCH_Q, ARCHerr_Q: the coefficients of the Q lagged squared
%       innovations, and their errors
% leverage, leverageerr: the leverage coefficient and its error (gjr only, otherwise NaN)
% distDoF: the degrees of freedom of a Student's t innovation distribution (NaN
%       otherwise)
% Goodness of fit:
% LLF, aic, bic: the log-likelihood, AIC and BIC per observation (divided by the
%       series length), of the maximum-likelihood fit
% summaryexitflag: 1 if the optimizer converged, otherwise 0
% persistence: the sum of the ARCH and GARCH coefficients (plus half the leverage
%       coefficient for 'gjr'), near 1 for near-integrated volatility
% uncondVar: the implied long-run (unconditional) variance (NaN if persistence is
%       0.999 or more)
% The fitted conditional variance (sigmas) through the series:
% maxsigma, minsigma, rangesigma, stdsigma, meansigma: its maximum, minimum, range,
%       standard deviation and mean
% Tests for remaining heteroskedasticity, comparing the (whitened) series and the
% residuals standardized by the conditional standard deviation, stde:
% engle_mean_diff_p, engle_max_diff_p: the mean and maximum, over lags 1 to 20, of
%       the change in p-value of Engle's ARCH test from the series to stde
% lbq_mean_diff_p, lbq_max_diff_p: the same for the Ljung-Box Q-test of the squared
%       series
% engle_pval_stde_1, engle_pval_stde_5, engle_pval_stde_10: the p-values of Engle's
%       ARCH test on stde at lags 1, 5 and 10
% minenglepval_stde, maxenglepval_stde: the minimum and maximum of those p-values
%       over lags 1 to 20
% lbq_pval_stde_1, lbq_pval_stde_5, lbq_pval_stde_10: the p-values of the Ljung-Box
%       Q-test on stde^2 at lags 1, 5 and 10
% minlbqpval_stde2, maxlbqpval_stde2: the minimum and maximum of those p-values over
%       lags 1 to 20
% ac1_stde2: the lag-1 autocorrelation of stde^2
% diff_ac1: the lag-1 autocorrelation of the squared series minus that of stde^2
% Summary of the standardized residuals, stde, from MF_ResidualAnalysis ('full'),
% with the prefix zres_:
% zres_meane, zres_meanabs, zres_stde, zres_maxonstd: mean, mean absolute value,
%       standard deviation, and largest absolute value (in standard deviations)
% zres_ac1, zres_ac2, zres_ac3: autocorrelation at lags 1 to 3
% zres_propbth: proportion of the autocorrelations at lags 1 to 25 within the
%       significance band +/- 2.6/sqrt(N)
% zres_ftbth: the first lag at which the autocorrelation is within that band
% zres_taurat: decorrelation time relative to that of the series
% zres_normksstat: Kolmogorov-Smirnov statistic against a normal distribution
% zres_sws, zres_swm: variation across 5 windows of the local standard deviation and
%       mean, relative to the overall standard deviation
% zres_popt, zres_minsbc: the order (1 to 10) of the AR model fitted to stde, chosen
%       by the Schwarz Bayesian criterion, and its criterion value
%
% ---NOTES:
% Conditional heteroskedasticity or not. The ARCH and GARCH coefficients are not
% identified when the series has no conditional heteroskedasticity: the likelihood is
% then flat in the GARCH coefficient, and the coefficients an optimizer reports depend on
% where it stops. The variance model is therefore reported only if it beats constant
% variance in a likelihood-ratio test; its statistic is rescaled by 2/(kurtosis - 1) of the
% standardized residuals for Gaussian innovations (heavy-tailed innovations inflate the
% unscaled statistic) and must reach 4 (6 for 'gjr', which has more room to fit noise). The
% null distribution of the statistic is not chi-squared, because a coefficient lies on the
% boundary of its range and another is not identified under the null (Andrews, Econometrica
% 69, 683 (2001), DOI: 10.1111/1468-0262.00210), so the critical values were calibrated by
% simulating white noise: 3% to 6% of white-noise series of 100 to 10000 samples exceed 4 (Gaussian,
% t(5), uniform and lognormal noise), and about 5% exceed 6 for 'gjr'. When the
% test is not passed, hasStructure is 0 and the constant-variance model is reported: the ARCH,
% GARCH and leverage coefficients and the persistence are 0, the conditional variance is the
% mean squared series, and the standard errors are NaN. The log-likelihood, AIC and BIC are
% always those of the fitted model, which is determined even where its coefficients are not.
%
% Student-t innovations. The degrees of freedom are reported only if the t fits better
% than a Gaussian with the same variance model in a likelihood-ratio test, whose null
% distribution is a 50:50 mixture of chi-squared(0) and chi-squared(1) because a
% Gaussian is the limit of infinite degrees of freedom, on the boundary (Self and Liang,
% J. Am. Stat. Assoc. 82, 605 (1987), DOI: 10.1080/01621459.1987.10478472): the statistic
% must reach 2.7055, the 5% point. Otherwise distDoF is 100, the upper bound of the fit
% (2.2 is the lower bound).
%
% Standard errors are from the outer product of the per-observation scores.
%
% Bounds. The constant of the variance equation is bounded below by 1e-6 (the variance of
% 0.1% of a standard deviation), and persistence above by 0.9999. For series with long
% stretches of (near) zero values (such as digital silence) the Gaussian likelihood has no
% maximum, because the variance can shrink without limit where the series is zero; the fit
% then sits at the bound (the constant within 1% of it). In that case LLF, aic and bic are NaN
% and, if the variance model is reported, so is everything computed from the conditional
% variance: the constant, the standard errors, uncondVar, the sigma summaries and the
% statistics of the standardized residuals. The coefficients, persistence, leverage and
% distDoF are still reported.
%
% Glosten, Jagannathan and Runkle's result for gjr: the leverage term acts on negative
% shocks only, so persistence = GARCH + ARCH + leverage/2. The leverage coefficient may be
% negative (down to -ARCH), so its sign is informative.
%
% Registered in hctsa's feature set: GARCH_1, ARCH_1, persistence, LLF, maxsigma, minsigma,
% engle_mean_diff_p, engle_max_diff_p, diff_ac1, engle_pval_stde_5,
% zres_meanabs, zres_maxonstd, zres_normksstat and zres_sws; and leverage for 'gjr' and
% distDoF for 't'. The other outputs are returned for inspection: the standard errors depend on
% the coefficients and the kurtosis rather than the structure and are invalid at a bound of the
% fit, and the rest repeat a registered output or describe the whitened series rather than the
% variance model (FeatureSets/renamed_operations.csv gives the reason for each).

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

beVocal = false; % Whether to display commentary on the fitting process

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(preproc)
	% Do an autoregressive preprocessing that maximizes stationarity/whitening
	preproc = 'ar';
end

% Fit what type of GARCH model?
if nargin < 3 || isempty(P)
	% Fit the default GARCH model
	P = 1;
end

if nargin < 4 || isempty(Q)
	% Fit the default GARCH model
	Q = 1;
end

% randomSeed: how to treat the randomization
if nargin < 5
	randomSeed = [];
end

if nargin < 6 || isempty(modelType)
	modelType = 'garch';
end

if nargin < 7 || isempty(innovationDist)
	innovationDist = 'gaussian';
end

if ~ismember(modelType, {'garch', 'gjr'})
	error('Unknown modelType ''%s'' (should be ''garch'' or ''gjr'')', modelType);
end
if ~ismember(innovationDist, {'gaussian', 't'})
	error('Unknown innovationDist ''%s'' (should be ''gaussian'' or ''t'')', innovationDist);
end
isGJR = strcmp(modelType, 'gjr');
isT = strcmp(innovationDist, 't');

% ------------------------------------------------------------------------------
%% (1) Data preprocessing
% ------------------------------------------------------------------------------
y = BF_Whiten(y, preproc, beVocal, randomSeed);

y = zscore(y); % z-score the time series (after whitening)

% Length of the (potentially whitened) time series, y
% Note that this could be different to the original (if choose a differencing, e.g.)
N = length(y);

% Now have the preprocessed time series saved over y.

% ------------------------------------------------------------------------------
%% (2) Data pre-estimation
% ------------------------------------------------------------------------------
% Statistics to compare to those of the residuals after modeling:
% (i) Engle's ARCH test, lags 1:20
Engle_pValue_y = BF_ARCHTest(y, 20);

% (ii) Ljung-Box Q-test of the squared series, lags 1:20
lbq_pValue_y2 = BF_LjungBoxTest(y.^2, 20);

% ------------------------------------------------------------------------------
%% (3) Fit the model
% ------------------------------------------------------------------------------
try
	D = BF_GARCHFit('data', y);
	gauss = BF_GARCHFit('mle', D, P, Q, isGJR, false);
	nullG = BF_GARCHFit('null', D, false);
	if isT
		sT = BF_GARCHFit('spec', D, P, Q, isGJR, true);
		alt = BF_GARCHFit('mle', D, P, Q, isGJR, true, ...
					BF_GARCHFit('tailstarts', gauss, sT), true);
		nullFit = BF_GARCHFit('null', D, true);
	else
		alt = gauss;
		nullFit = nullG;
	end
	if ~isfinite(alt.llf) || any(~isfinite(alt.theta))
		error('non-finite estimates');
	end
	% Does the variance model beat constant variance?
	[~, hasStructure] = BF_GARCHFit('structure', alt, nullFit, D, isGJR);
catch emsg
	error('GARCH fit failed (data does not allow a valid GARCH model to be estimated): %s', emsg.message);
end

% If not, the variance model is not identified: report the constant-variance model
if hasStructure
	fit = alt;
else
	fit = nullFit;
end
S = fit.spec;
theta = fit.theta;

% ------------------------------------------------------------------------------
%% (4) Return statistics on fit
% ------------------------------------------------------------------------------
out.hasStructure = double(hasStructure);
out.offset = 0; % (zero-mean model of a z-scored series)
out.constant = theta(1);

% Coefficients and their standard errors
if hasStructure
	errors = BF_GARCHFit('errors', theta, D, S);
	out.constanterr = errors(1);
	for i = 1:P
		out.(sprintf('GARCH_%u', i)) = theta(S.ib(i));
		out.(sprintf('GARCHerr_%u', i)) = errors(S.ib(i));
	end
	for i = 1:Q
		out.(sprintf('ARCH_%u', i)) = theta(S.ia(i));
		out.(sprintf('ARCHerr_%u', i)) = errors(S.ia(i));
	end
	if isGJR
		out.leverage = theta(S.ig(1));
		out.leverageerr = errors(S.ig(1));
		persistence = sum(theta(S.ib)) + sum(theta(S.ia)) + sum(theta(S.ig))/2;
	else
		out.leverage = NaN;
		out.leverageerr = NaN;
		persistence = sum(theta(S.ib)) + sum(theta(S.ia));
	end
else
	out.constanterr = NaN;
	for i = 1:P
		out.(sprintf('GARCH_%u', i)) = 0;
		out.(sprintf('GARCHerr_%u', i)) = NaN;
	end
	for i = 1:Q
		out.(sprintf('ARCH_%u', i)) = 0;
		out.(sprintf('ARCHerr_%u', i)) = NaN;
	end
	if isGJR
		out.leverage = 0;
	else
		out.leverage = NaN;
	end
	out.leverageerr = NaN;
	persistence = 0;
end

% Degrees of freedom of the Student's t distribution
if isT
	% Student-t against Gaussian innovations with the same variance model
	if hasStructure
		lrNu = 2 * (alt.llf - gauss.llf);
	else
		lrNu = 2 * (nullFit.llf - nullG.llf);
	end
	if lrNu >= BF_GARCHFit('taillr')
		out.distDoF = 1 / theta(S.ik);
	else
		out.distDoF = 100; % the upper bound: not distinguishable from Gaussian
	end
else
	out.distDoF = NaN;
end

% More statistics given from the fit
% (the likelihood of the fitted model is determined even where its parameters are not)
out.LLF = alt.llf / N; % log-likelihood per observation (the total scales with length)
out.summaryexitflag = double(alt.dec < 1e-8); % whether the optimizer converged
out.aic = (-2*alt.llf + 2*alt.k) / N; % per observation, as for LLF
out.bic = (-2*alt.llf + alt.k*log(N)) / N;

% Persistence (sum of ARCH + GARCH coefficients, i.e. how long volatility shocks
% persist) and the implied long-run (unconditional) variance. For gjr, the leverage term
% acts on negative shocks only, so on average half the time (Glosten, Jagannathan and
% Runkle, 1993): persistence = GARCH + ARCH + leverage/2.
out.persistence = persistence;
if persistence < 0.999
	out.uncondVar = theta(1) / (1 - persistence);
else
	out.uncondVar = NaN;
end

% ------------------------------------------------------------------------------
%% Sigmas, the time series of conditional variances
% ------------------------------------------------------------------------------
% For a time series with strong ARCH/GARCH effects, this will fluctuate;
% otherwise will be flat
sigmas = fit.h;
out.maxsigma = max(sigmas);
out.minsigma = min(sigmas);
out.rangesigma = max(sigmas) - min(sigmas); % very similar information to max(sigma) for most time series
out.stdsigma = std(sigmas);
out.meansigma = mean(sigmas);
if ~hasStructure
	out.rangesigma = 0; % (exactly zero for constant variance)
	out.stdsigma = 0;
end

% ------------------------------------------------------------------------------
%% Check residuals
% ------------------------------------------------------------------------------
res = (0 - y); % residuals (mean process minus data, the MF_ResidualAnalysis convention)
stde = res ./ sqrt(sigmas); % standardize residuals by conditional standard deviation
stde2 = stde.^2;

Engle_pValue_stde = BF_ARCHTest(stde, 20);
lbq_pValue_stde2 = BF_LjungBoxTest(stde2, 20);

% Ok, so now we've corrected for GARCH effects, how does this 'improve' the
% randomness/correlation in our signal. If the signal is much less
% correlated now, it is a signature that GARCH effects were significant in
% the original signal

% Mean/max improvement in Engle pValue
out.engle_mean_diff_p = mean(Engle_pValue_stde - Engle_pValue_y);
out.engle_max_diff_p = max(Engle_pValue_stde - Engle_pValue_y);

% Mean/max improvement in lbq pValue for squared time series
out.lbq_mean_diff_p = mean(lbq_pValue_stde2 - lbq_pValue_y2);
out.lbq_max_diff_p = max(lbq_pValue_stde2 - lbq_pValue_y2);

% Raw values:
out.engle_pval_stde_1 = Engle_pValue_stde(1);
out.engle_pval_stde_5 = Engle_pValue_stde(5);
out.engle_pval_stde_10 = Engle_pValue_stde(10);
out.minenglepval_stde = min(Engle_pValue_stde);
out.maxenglepval_stde = max(Engle_pValue_stde);

out.lbq_pval_stde_1 = lbq_pValue_stde2(1);
out.lbq_pval_stde_5 = lbq_pValue_stde2(5);
out.lbq_pval_stde_10 = lbq_pValue_stde2(10);
out.minlbqpval_stde2 = min(lbq_pValue_stde2);
out.maxlbqpval_stde2 = max(lbq_pValue_stde2);

% Use MF_ResidualAnalysis on the standardized innovations
% 1) Get statistics on the standardized innovations, prefixed zres_
residout = MF_ResidualAnalysis(stde, y, 'full');

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields)
	out.(sprintf('zres_%s', fields{k})) = residout.(fields{k});
end

out.ac1_stde2 = CO_AutoCorr(stde2, 1, 'Fourier');
out.diff_ac1 = CO_AutoCorr(y.^2, 1, 'Fourier') - CO_AutoCorr(stde2, 1, 'Fourier');

% ------------------------------------------------------------------------------
%% A fit with the variance constant at its floor has no likelihood maximum
% ------------------------------------------------------------------------------
% What depends on the likelihood or on the variance path is then not a property of the series
if BF_GARCHFit('atfloor', alt.theta)
	out.LLF = NaN;
	out.aic = NaN;
	out.bic = NaN;
	if hasStructure
		keepFields = {'hasStructure', 'offset', 'leverage', 'distDoF', 'persistence', 'summaryexitflag'};
		fields = fieldnames(out);
		for k = 1:length(fields)
			if ~ismember(fields{k}, keepFields) && isempty(regexp(fields{k}, '^(GARCH|ARCH)_', 'once'))
				out.(fields{k}) = NaN;
			end
		end
	end
end

end
