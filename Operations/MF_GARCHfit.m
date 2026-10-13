function out = MF_GARCHfit(y, preproc, P, Q, randomSeed, modelType, innovationDist)
% MF_GARCHfit   A GARCH model of the changing variance of the series, and what it leaves in the residuals.
%
% Simulates a procedure for fitting Generalized Autoregressive Conditional
% Heteroskedasticity (GARCH) models to a time series, namely:
%
% (1) Preprocessing the data to remove strong trends (and, optionally, to whiten it),
% (2) Pre-estimation to calculate initial correlation properties of the time
%       series and motivate a GARCH model,
% (3) Fitting a GARCH model, returning goodness-of-fit statistics and parameters
%           of the fitted model, and
% (4) Post-estimation, involving statistics on the residuals and the standardized
%           residuals.
%
% The preprocessing is BF_Whiten: the series is detrended and, for preproc = 'ar',
% replaced by the preprocessing (from PP_PreProcess) that maximizes whiteness under an
% AR(2) model, if one improves on the original series by more than 5%. The result is
% then z-scored, so the model is of the variance around a zero mean.
%
% Uses functions from MATLAB's Econometrics Toolbox: garch, gjr, egarch, estimate,
% infer, archtest, lbqtest, autocorr, parcorr, aicbic.
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
%               symmetric response to shocks), 'gjr' (GJR-GARCH, adds a
%               leverage/asymmetry term so negative and positive shocks can
%               have different effects on variance), or 'egarch'
%               (exponential GARCH, models log-variance, also asymmetric).
%
% innovationDist, the assumed innovation distribution: 'gaussian' (default)
%               or 't' (Student's t, estimates a degrees-of-freedom
%               parameter to capture fat tails beyond what GARCH-filtering
%               alone accounts for).
%
% ---OUTPUTS:
% Fitted parameters (standard errors have 'err' in the name):
% constant, constanterr: the constant term of the variance equation, and its error
% offset: the mean offset of the model (0 for a z-scored series)
% GARCH_1, GARCHerr_1, ..., GARCH_P, GARCHerr_P: the coefficients of the P lagged
%       variances, and their errors (NaN if that lag was dropped from the fit)
% ARCH_1, ARCHerr_1, ..., ARCH_Q, ARCHerr_Q: the coefficients of the Q lagged squared
%       innovations, and their errors
% leverage, leverageerr: the leverage coefficient and its error (gjr/egarch only,
%       otherwise NaN)
% distDoF: the degrees of freedom of a Student's t innovation distribution (NaN
%       otherwise)
% Goodness of fit:
% LLF, aic, bic: the log-likelihood, AIC and BIC per observation (divided by the
%       series length)
% summaryexitflag: the exit flag of the fitting procedure
% persistence: the sum of the ARCH and GARCH coefficients (plus half the leverage
%       coefficient for 'gjr'; NaN for 'egarch'), near 1 for near-integrated volatility
% uncondVar: the implied long-run (unconditional) variance (NaN if persistence is
%       0.999 or more, or for 'egarch')
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
% Only the P = 1, Q = 1 registration remains. The P = 1, Q = 2 registration was
% dropped as almost entirely redundant with it (r = 0.9-1.0 across nearly every
% output field, on two collections of real-world series). MF_GARCHcompare, which
% varies P and Q over a grid, is unaffected.
%
% persistence and uncondVar are standard GARCH diagnostics (persistence near 1
% signals near-integrated, IGARCH-like volatility clustering). uncondVar is set to
% NaN near the boundary because the model's value is numerically meaningless there.
% For 'gjr', persistence includes half the leverage coefficient (Glosten,
% Jagannathan and Runkle's result: the leverage term acts on negative shocks only).
%
% The modelType and innovationDist arguments, and the leverage, leverageerr and
% distDoF outputs, cover asymmetric volatility response and fat-tailed innovations.
% Registered variants: MF_GARCHfit_ar_P1_Q1_gjr and MF_GARCHfit_ar_P1_Q1_t each
% change one thing from the P1_Q1 baseline. 'egarch' is supported but not registered
% (its ARCH coefficient hit an apparent boundary of 1.0 in 2 of 3 real fits,
% unexplained).

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

% Check that an Econometrics Toolbox license is available:
BF_CheckToolbox('econometrics_toolbox');

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

% ------------------------------------------------------------------------------
%% (1) Data preprocessing
% ------------------------------------------------------------------------------
% Save the original, unprocessed time series
y0 = y;

y = BF_Whiten(y, preproc, beVocal, randomSeed);

y = zscore(y); % z-score the time series (after whitening)

% Length of the (potentially whitened) time series, y
% Note that this could be different to the original, y0 (if choose a differencing, e.g.)
N = length(y);

% Now have the preprocessed time series saved over y.
% The original, unprocessed time series is retained in y0.
% (Note that y=y0 is possible; when all preprocessings are found to be
%   worse at the given criterion).

% ------------------------------------------------------------------------------
%% (2) Data pre-estimation
% ------------------------------------------------------------------------------
% Aim is to return some statistics indicating the suitability of this class
% of modeling.
% Will use the statistics to motivate a GARCH model in the next
% section.
% Will use the statistics to compare to features of the residuals after
% modeling.

% (i) Engle's ARCH test
%       look at autoregressive lags 1:20
%       use the 10% significance level
[Engle_h_y, Engle_pValue_y, Engle_stat_y, Engle_cValue_y] = archtest(y, 'lags', 1:20, 'alpha', 0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2, 'lags', 1:20, 'alpha', 0.1);

% (Commented out as a historical marker: unused by any output; saves wasted compute,
%  as in MF_GARCHcompare.)
% (iii) Correlation in time series: autocorrelation
% [ACF_y, Lags_acf_y, bounds_acf_y] = autocorr(y, 'NumLags', 20);
% [ACF_var_y, Lags_acf_var_y, bounds_acf_var_y] = autocorr(y.^2, 'NumLags', 20);

% (iv) Partial autocorrelation function: PACF
% [PACF_y, Lags_pacf_y, bounds_pacf_y] = parcorr(y, 'NumLags', 20);

% ------------------------------------------------------------------------------
%% (3) Create an appropriate GARCH model
% ------------------------------------------------------------------------------
switch modelType
case 'garch'
	GModel = garch(P, Q); % GARCH degree P, ARCH degree Q
case 'gjr'
	GModel = gjr(P, Q); % adds a leverage/asymmetry term
case 'egarch'
	GModel = egarch(P, Q); % log-variance, also asymmetric
otherwise
	error('Unknown modelType ''%s'' (should be ''garch'', ''gjr'', or ''egarch'')', modelType);
end

% Include a constant in the GARCH model
GModel.Constant = NaN;

if ~strcmp(innovationDist, 'gaussian')
	GModel.Distribution = innovationDist; % e.g., 't' estimates a degrees-of-freedom parameter
end

% Fit the model
try
	[Gfit, estParamCov, LLF, info] = estimate(GModel, y, 'Display', 'off');
	% Estimate standard errors using variance/covariance matrix:
	errors = sqrt(diag(estParamCov));
	% [coeff, errors, LLF, innovations, sigmas, summary] = garchfit(GModel,y);
catch emsg
	error('GARCH fit failed (data does not allow a valid GARCH model to be estimated): %s', emsg.message);
	% Sometimes this happens for some time series (e.g., when it removes some GARCH
	% lags and makes the resulting model invalid)
end

% ------------------------------------------------------------------------------
%% (4) Return statistics on fit
% ------------------------------------------------------------------------------

% (i) Return coefficients, and their standard errors as seperate statistics
% ___Mean_Process___
% --Constant--
if isprop(Gfit, 'Constant')
	out.constant = Gfit.Constant;
	out.constanterr = errors(1);
end

% __Variance_Process___
% -- Offset (should be zero for z-scored time series)--
if isprop(Gfit, 'Offset')
	out.offset = Gfit.Offset;
end

% The variance-covariance matrix from estimate has one row/column per estimated parameter
% (constant, GARCH lags, ARCH lags, then any leverage/DoF), in that order, whether or
% not a coefficient was estimated at exactly zero, so the error for a lag is simply
% indexed by its position.

% -- GARCH component --
for i = 1:P
	if isprop(Gfit, 'GARCH') && length(Gfit.GARCH) >= i
		out.(sprintf('GARCH_%u', i)) = Gfit.GARCH{i};
		% New (in this way shit) format means that this no longer works for
		% custom GARCH models (you can no longer index a particular
		% error) ///
		if Gfit.GARCH{i} == 0
			% no fit at this lag, even though it was specified
			out.(sprintf('GARCHerr_%u', i)) = NaN; % first is the constant
		else
			out.(sprintf('GARCHerr_%u', i)) = errors(1 + i); % first is the constant
		end
	else
		% fitted GARCH model not as specified
		out.(sprintf('GARCH_%u', i)) = NaN;
		out.(sprintf('GARCHerr_%u', i)) = NaN; % first is the constant
	end
end

% -- ARCH component --
for i = 1:Q
	if isprop(Gfit, 'ARCH') && length(Gfit.ARCH) >= i
		out.(sprintf('ARCH_%u', i)) = Gfit.ARCH{i};
		if Gfit.ARCH{i} == 0
			% No fit at this specified lag
			out.(sprintf('ARCHerr_%u', i)) = NaN; % constant, then GARCH, then ARCH
		else
			out.(sprintf('ARCHerr_%u', i)) = errors(1 + length(Gfit.GARCH) + i); % constant, then GARCH, then ARCH
		end
	else
		% ARCH fit not as specified
		out.(sprintf('ARCH_%u', i)) = NaN;
		out.(sprintf('ARCHerr_%u', i)) = NaN;
	end
end

% -- Leverage/asymmetry component (gjr/egarch only) --
% For the single-lag models this operation registers, the extra
% parameter (leverage or DoF, below) is always the LAST element of `errors`.
if isprop(Gfit, 'Leverage') && ~isempty(Gfit.Leverage)
	out.leverage = Gfit.Leverage{1};
	out.leverageerr = errors(end);
else
	out.leverage = NaN;
	out.leverageerr = NaN;
end

% -- Innovation-distribution degrees of freedom (Student's t only) --
if strcmp(Gfit.Distribution.Name, 't')
	out.distDoF = Gfit.Distribution.DoF;
else
	out.distDoF = NaN;
end

% More statistics given from the fit
out.LLF = LLF / N; % log-likelihood per observation (the total scales with length)

out.summaryexitflag = info.exitflag; % whether the fit worked ok.
% This is just a record, really, since the numerical values are only
% symbolic.

nparams = sum(any(estParamCov)); % number of parameters
% (Not registered as a feature: for a fixed P/Q/AR order, this is a structural property of
%  the model specification, not of the data, so it is constant across series -- confirmed
%  on both validation datasets. Still computed here for the AIC/BIC below.)

% use aicbic function
[AIC, BIC] = aicbic(LLF, nparams, N); % aic and bic of fit
out.aic = AIC / N; % per observation, as for LLF
out.bic = BIC / N;

% Persistence (sum of ARCH + GARCH coefficients, i.e. how long volatility
% shocks persist) and implied long-run (unconditional) variance. Persistence
% close to 1 indicates near-integrated (IGARCH-like) volatility clustering;
% >= 1 would mean no finite unconditional variance exists. In practice the
% estimate() optimizer enforces a stationarity constraint with a small
% internal tolerance, so near-boundary fits land just under 1 (observed
% exactly 0.9999998 on real data) rather than at/over it -- uncondVar is
% NaN'd out here rather than only guarding on persistence>=1, since a
% denominator that small makes the value numerically meaningless (dominated
% by the optimizer's boundary tolerance, not the data) well before
% persistence formally reaches 1.
%
% For asymmetric (gjr) models, persistence isn't just GARCH+ARCH: the
% leverage term only applies on negative shocks, so under the fitted
% (symmetric, zero-mean) innovation distribution it contributes on average
% half the time -- persistence = GARCH+ARCH+Leverage/2 (Glosten-Jagannathan-
% Runkle 1993). Confirmed empirically: this matches the model object's own
% UnconditionalVariance property almost exactly (to ~4 sig figs) across 5
% real series, whereas omitting the Leverage/2 term gives nonsense,
% including NEGATIVE "variances" on 2/5 series. uncondVar itself is read
% directly from Gfit.UnconditionalVariance (native to garch/gjr/egarch
% objects) rather than hand-derived, to avoid this class of formula bug --
% still requires the persistence-based NaN guard above near the boundary,
% since the native property blows up there too, not just a hand-rolled one.
% egarch's log-variance recursion doesn't reduce to a simple coefficient
% sum, so persistence/uncondVar are left NaN there (egarch is unregistered
% currently anyway).
switch modelType
case 'garch'
	out.persistence = sum(cellfun(@(c) c, Gfit.GARCH)) + sum(cellfun(@(c) c, Gfit.ARCH));
case 'gjr'
	out.persistence = sum(cellfun(@(c) c, Gfit.GARCH)) + sum(cellfun(@(c) c, Gfit.ARCH)) + ...
						sum(cellfun(@(c) c, Gfit.Leverage))/2;
otherwise % egarch
	out.persistence = NaN;
end
if ~isnan(out.persistence) && out.persistence < 0.999
	out.uncondVar = Gfit.UnconditionalVariance;
else
	out.uncondVar = NaN;
end

% ------------------------------------------------------------------------------
%% Sigmas, the time series of conditional variances
% ------------------------------------------------------------------------------
% Estimate it:
[sigmas, logL] = infer(Gfit, y);
% For a time series with strong ARCH/GARCH effects, this will fluctuate;
% otherwise will be quite flat...
out.maxsigma = max(sigmas);
out.minsigma = min(sigmas);
out.rangesigma = max(sigmas) - min(sigmas); % very similar information to max(sigma) for most time series
out.stdsigma = std(sigmas);
out.meansigma = mean(sigmas);

% ------------------------------------------------------------------------------
%% Check residuals
% ------------------------------------------------------------------------------
res = (Gfit.Offset - y); % residuals (mean process minus data, the MF_ResidualAnalysis convention)
stde = res ./ sqrt(sigmas); % standardize residuals by conditional standard deviation
stde2 = stde.^2;

% (i) Engle's ARCH test
%       look at autoregressive lags 1:20
%       use the 10% significance level
[Engle_h_stde, Engle_pValue_stde, Engle_stat_stde, Engle_cValue_stde] = archtest(stde, 'lags', 1:20, 'alpha', 0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_stde2, lbq_pValue_stde2, lbq_stat_stde2, lbq_cValue_stde2] = lbqtest(stde2, 'lags', 1:20, 'alpha', 0.1);

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

% (iii) Correlation in time series: autocorrelation
% autocorrs_y = CO_AutoCorr(y,1:20);
% autocorrs_var = CO_AutoCorr(y.^2,1:20);
% [ACF_y,Lags_acf_y,bounds_acf_y] = autocorr(e,20,[],[]);
% [ACF_var,Lags_acf_var,bounds_acf_var] = autocorr(e.^2,20,[],[]);

% (iv) Partial autocorrelation function: PACF
% [PACF,Lags_pacf,bounds_pacf] = parcorr(e,20,[],[]);

% Use MF_ResidualAnalysis on the standardized innovations
% 1) Get statistics on the standardized innovations, prefixed zres_ (the old stde_
%    prefix produced field names like stde_stde)
residout = MF_ResidualAnalysis(stde, y, 'full');

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields);
	out.(sprintf('zres_%s', fields{k})) = residout.(fields{k});
end

out.ac1_stde2 = CO_AutoCorr(stde2, 1, 'Fourier');
out.diff_ac1 = CO_AutoCorr(y.^2, 1, 'Fourier') - CO_AutoCorr(stde2, 1, 'Fourier');

%% (5) Comparison to other models
% e.g., does the additional heteroskedastic component improve the model fit
% over just the conditional mean component of the model.

end
