function out = MF_GARCHfit(y,preproc,P,Q,randomSeed)
% MF_GARCHfit   GARCH time-series modeling.
%
% Simulates a procedure for fitting Generalized Autoregressive Conditional
% Heteroskedasticity (GARCH) models to a time series, namely:
%
% (1) Preprocessing the data to remove strong trends,
% (2) Pre-estimation to calculate initial correlation properties of the time
%       series and motivate a GARCH model,
% (3) Fitting a GARCH model, returning goodness of fit statistics and parameters
%           of the fitted model, and
% (4) Post-estimation, involves calculating statistics on residuals and
%           standardized residuals.
%
% The idea is that all of these stages can be pre-specified or skipped using
% arguments to the function.
%
% Uses functions from Matlab's Econometrics Toolbox: archtest, lbqtest,
% autocorr, parcorr, garchset, garchfit, garchcount, aicbic
%
% All methods implemented are from Matlab's Econometrics Toolbox, including
% Engle's ARCH test (archtest), the Ljung-Box Q-test (lbqtest), estimating the
% partial autocorrelation function (parcorr), as well as specifying (garchset)
% and fitting (garchfit) the GARCH model to the time series.
%
% As part of this code, a very basic automatic pre-processing routine,
% PP_ModelFit, is implemented, that applies a range of pre-processings and
% returns the preprocessing of the time series that shows the worst fit to an
% AR(2) model.
%
% In the case that no simple transformations of the time series are
% significantly more stationary/less trivially correlated than the original time
% series (by more than 5%), the original time series is simply used without
% transformation.
%
% Where r and m are the autoregressive and moving average orders of the model,
% respectively, and p and q control the conditional variance parameters.
%
%---INPUTS:
% y, the input time series
%
% preproc, the preprocessing to apply, can be 'ar' or 'none'
%
% P, the GARCH model order
%
% Q, the ARCH model order
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (for pre-processing: PP_PreProcess)

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

beVocal = 0; % Whether to display commentary on the fitting process

% Check that an Econometrics Toolbox license is available:
BF_CheckToolbox('econometrics_toolbox')

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

% ------------------------------------------------------------------------------
%% (1) Data preprocessing
% ------------------------------------------------------------------------------
% Save the original, unprocessed time series
y0 = y;

y = BF_Whiten(y,preproc,beVocal,randomSeed);

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
[Engle_h_y, Engle_pValue_y, Engle_stat_y, Engle_cValue_y] = archtest(y,'lags',1:20,'alpha',0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,'lags',1:20,'alpha',0.1);


% (iii) Correlation in time series: autocorrelation
[ACF_y, Lags_acf_y, bounds_acf_y] = autocorr(y,20,[],[]);
[ACF_var_y, Lags_acf_var_y, bounds_acf_var_y] = autocorr(y.^2,20,[],[]);

% (iv) Partial autocorrelation function: PACF
[PACF_y, Lags_pacf_y, bounds_pacf_y] = parcorr(y,20,[],[]);


% ------------------------------------------------------------------------------
%% (3) Create an appropriate GARCH model
% ------------------------------------------------------------------------------
GModel = garch(P,Q); % ARCH order P, GARCH order Q

% Include a constant in the GARCH model
GModel.Constant = NaN;

% Fit the model
try
    [Gfit, estParamCov, LLF, info] = estimate(GModel,y,'Display','off');
    % Estimate standard errors using variance/covariance matrix:
    errors = sqrt(diag(estParamCov));
    % [coeff, errors, LLF, innovations, sigmas, summary] = garchfit(GModel,y);
catch emsg
    error('GARCH fit failed (data does not allow a valid GARCH model to be estimated): %s',emsg.message);
    % Sometimes this happens for some time series (e.g., when it removes some GARCH
    % lags and makes the resulting model invalid)
end

% ------------------------------------------------------------------------------
%% (4) Return statistics on fit
% ------------------------------------------------------------------------------

% (i) Return coefficients, and their standard errors as seperate statistics
% ___Mean_Process___
% --Constant--
if isprop(Gfit,'Constant')
    out.constant = Gfit.Constant;
    out.constanterr = errors(1);
end

% __Variance_Process___
% -- Offset (should be zero for z-scored time series)--
if isprop(Gfit,'Offset')
    out.offset = Gfit.Offset;
end


indexAdjust = 0; % required because sometimes you fit at fewer lags than you
                 % specified, but the errors output is a vector,
                 % so sadly you have to keep count...

% -- GARCH component --
for i = 1:P
    if isprop(Gfit,'GARCH') && length(Gfit.GARCH)>=i
        out.(sprintf('GARCH_%u',i)) = Gfit.GARCH{i};
        % New (in this way shit) format means that this no longer works for
        % custom GARCH models (you can no longer index a particular
        % error) ///
        if Gfit.GARCH{i}==0
            % no fit at this lag, even though it was specified
            indexAdjust = indexAdjust + 1;
            out.(sprintf('GARCHerr_%u',i)) = NaN; % first is the constant
        else
            out.(sprintf('GARCHerr_%u',i)) = errors(1+i-indexAdjust); % first is the constant
        end
    else
        % fitted GARCH model not as specified
        out.(sprintf('GARCH_%u',i)) = NaN;
        out.(sprintf('GARCHerr_%u',i)) = NaN; % first is the constant
    end
end

% -- ARCH component --
for i = 1:Q
    if isprop(Gfit,'ARCH') && length(Gfit.ARCH)>=i
        out.(sprintf('ARCH_%u',i)) = Gfit.ARCH{i};
        if Gfit.ARCH{i}==0
            % No fit at this specified lag
            out.(sprintf('ARCHerr_%u',i)) = NaN; % constant, then GARCH, then ARCH
            indexAdjust = indexAdjust + 1;
        else
            out.(sprintf('ARCHerr_%u',i)) = errors(1+length(Gfit.GARCH)+i-indexAdjust); % constant, then GARCH, then ARCH
        end
    else
        % ARCH fit not as specified
        out.(sprintf('ARCH_%u',i)) = NaN;
        out.(sprintf('ARCHerr_%u',i)) = NaN;
    end
end

% More statistics given from the fit
out.LLF = LLF; % log-likelihood function

out.summaryexitflag = info.exitflag; % whether the fit worked ok.
% This is just a record, really, since the numerical values are only
% symbolic.

nparams = sum(any(estParamCov)); % number of parameters
out.nparams = nparams;

% use aicbic function
[AIC, BIC] = aicbic(LLF,nparams,N); % aic and bic of fit
out.aic = AIC;
out.bic = BIC;

% ------------------------------------------------------------------------------
%% Sigmas, the time series of conditional variances
% ------------------------------------------------------------------------------
% Estimate it:
[sigmas,logL] = infer(Gfit,y);
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
res = (y - Gfit.Offset); % residuals (departures from mean process)
stde = res./sqrt(sigmas); % standardize residuals by conditional standard deviation
stde2 = stde.^2;

% (i) Engle's ARCH test
%       look at autoregressive lags 1:20
%       use the 10% significance level
[Engle_h_stde, Engle_pValue_stde, Engle_stat_stde, Engle_cValue_stde] = archtest(stde,'lags',1:20,'alpha',0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_stde2, lbq_pValue_stde2, lbq_stat_stde2, lbq_cValue_stde2] = lbqtest(stde2,'lags',1:20,'alpha',0.1);


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
% 1) Get statistics on standardized innovations
residout = MF_ResidualAnalysis(stde);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k = 1:length(fields);
    out.(sprintf('stde_%s',fields{k})) = residout.(fields{k});
end

out.ac1_stde2 = CO_AutoCorr(stde2,1,'Fourier');
out.diff_ac1 = CO_AutoCorr(y.^2,1,'Fourier') - CO_AutoCorr(stde2,1,'Fourier');


%% (5) Comparison to other models
% e.g., does the additional heteroskedastic component improve the model fit
% over just the conditional mean component of the model.


end
