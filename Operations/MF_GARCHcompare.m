function out = MF_GARCHcompare(y,preProc,pr,qr,randomSeed,beVocal)
% MF_GARCHcompare   Comparison of GARCH time-series models
%
% This code fits a set of GARCH(p,q) models to the time series and
% returns statistics on the goodness of fits across a range of p and
% q parameters.
%
% Uses the following functions from Matlab's Econometrics Toolbox: archtest,
% lbqtest, autocorr, parcorr, garchset, garchfit, garchcount, aicbic
%
% Compares all combinations of p and q and output statistics are on the models
% with the best fit.
%
% This operation focuses on the GARCH/variance component, and therefore
% attempts to pre-whiten and assumes a constant mean process (applies a linear
% detrending).
%
%---INPUTS:
% y, the input time series
% preProc, a preprocessing to apply:
%           (i) 'none': no preprocessing is performed
%           (ii) 'ar': performs a preprocessing that maximizes AR(2) whiteness,
%
% pr, a vector of model orders, p, to compare
%
% qr, a vector of model orders, q, to compare
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%
%---OUTPUTS: include log-likelihoods, Bayesian Information  Criteria (BIC),
% Akaike's Information Criteria (AIC), outputs from Engle's ARCH test and the
% Ljung-Box Q-test, and estimates of optimal model orders.

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
%% Check that an Econometrics Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('econometrics_toolbox')

% ------------------------------------------------------------------------------
%% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(preProc)
    preProc = 'none';
end

% GARCH (p) ARCH (q) parameters
if nargin < 3 || isempty(pr)
    pr = (1:3); % i.e., GARCH(1:3,qr)
end

if nargin < 4 || isempty(qr)
    qr = (1:3); % i.e., GARCH(pr,1:3);
end

% randomSeed: how to treat the randomization (in BF_Whiten)
if nargin < 5
    randomSeed = [];
end

% beVocal: whether to speak to the command line
if nargin < 6
    beVocal = 0; % (no by default)
end

% ------------------------------------------------------------------------------
%% (1) Data preprocessing
% ------------------------------------------------------------------------------
y0 = y; % the original, unprocessed time series

y = BF_Whiten(y,preProc,0,randomSeed); % Use a basic preprocessing to whiten the time series

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
y = zscore(y); % make sure the (whitened) time series is z-scored
N = length(y); % could be different to original (e.g., if chose a differencing above)

% Now have the preprocessed time series saved over y.
% The original, unprocessed time series is retained in y0.
% (Note that y=y0 is possible; when all preprocessings are found to be
%   worse at the given criterion).

%% (2) Data pre-estimation
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
% [Engle_h_y, Engle_pValue_y, Engle_stat_y, Engle_cValue_y] = archtest(y,1:20,0.1); % depricated syntax

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,'lags',1:20,'alpha',0.1);
% [lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,1:20,0.1); % depricated syntax
% [lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,1:20,0.1,[]); % depricated syntax


% (iii) Correlation in time series: autocorrelation
% autocorrs_y = CO_AutoCorr(y,1:20);
% autocorrs_var = CO_AutoCorr(y.^2,1:20);
[ACF_y, Lags_acf_y, bounds_acf_y] = autocorr(y,20,[],[]);
[ACF_var_y, Lags_acf_var_y, bounds_acf_var_y] = autocorr(y.^2,20,[],[]);

% (iv) Partial autocorrelation function: PACF
[PACF_y, Lags_pacf_y, bounds_pacf_y] = parcorr(y,20,[],[]);

% ------------------------------------------------------------------------------
%% (3) Create an appropriate GARCH model
% ------------------------------------------------------------------------------
% initialize statistics
np = length(pr);
nq = length(qr);
LLFs = NaN*ones(np,nq); % log-likelihood
AICs = NaN*ones(np,nq); % AIC
BICs = NaN*ones(np,nq); % BIC
Ks = NaN*ones(np,nq); % constant term in variance
meanarchps = NaN*ones(np,nq); % mean over 20 lags from Engle's ARCH test on standardized innovations
maxarchps = NaN*ones(np,nq); % maximum p-value over 20 lags from Engle's ARCH test on standardized innovations
meanlbqps = NaN*ones(np,nq); % mean lbq p-value over 20 lags from Q-test on squared standardized innovations
maxlbqps = NaN*ones(np,nq); % maximum p-value over 20 lags from Q-test on squared standardized innovations
isBad = zeros(np,nq);

for i = 1:np
    p = pr(i); % garch order
    for j = 1:nq
       q = qr(j); % arch order

       % (i) specify a zero-mean, Gaussian innovation GARCH(P,Q) model.
       GModel = garch(p,q);
       GModel.Constant = NaN;

       % (ii) fit the model
       [Gfit, estParamCov, LLF, info] = estimate(GModel,y,'Display','off');
       % [coeff, errors, LLF, innovations, sigmas, summary] = garchfit(spec,y);

       numParams = sum(any(estParamCov)); % number of parameters
       if numParams < p + q + 1
           if beVocal
               fprintf(1,'Bad fit at p = %u, q = %u\n',p,q);
           end
           isBad(i,j) = 1;
           continue; % didn't fit successfully; everything stays NaN
       end

       % (iii) store derived statistics on the fitted model
       LLFs(i,j) = LLF;
       [AIC, BIC] = aicbic(LLF,numParams,N); % aic and bic of fit
       AICs(i,j) = AIC;
       BICs(i,j) = BIC;
       Ks(i,j) = Gfit.Constant;

       % (iv) store summaries of standardized errors
       %    (less need to compare to original because now we're only
       %    comparing *between* models, which should all have the same
       %    baseline...?) It's just that the absolute values of these
       %    summaries are less meaningful; should be with respect to the
       %    original time series
       innovations = (y - Gfit.Offset); % residuals (departures from mean process)
       [sigmas,logL] = infer(Gfit,y); % estimate time series of conditional variance, sigmas
       stde = innovations./sqrt(sigmas); % standardized residuals
       stde2 = stde.^2; % squared

       % (i) Engle's ARCH test
       %       look at autoregressive lags 1:20
       [Engle_h_stde, Engle_pValue_stde, Engle_stat_stde, Engle_cValue_stde] = archtest(stde,'lags',1:20,'alpha',0.1);

       meanarchps(i,j) = mean(Engle_pValue_stde);
       maxarchps(i,j) = max(Engle_pValue_stde);

       % (ii) Ljung-Box Q-test
       %       look at autocorrelation at lags 1:20
       %       use the 10% significance level
       %       departure from randomness hypothesis test
       [lbq_h_stde2, lbq_pValue_stde2, lbq_stat_stde2, lbq_cValue_stde2] = lbqtest(stde2,'lags',1:20,'alpha',0.1);

       meanlbqps(i,j) = mean(lbq_pValue_stde2);
       maxlbqps(i,j) = max(lbq_pValue_stde2);

       % Difficult to take parameter values since their number if
       % changing...
    end
end

% Check if fits went well...
if all(isBad(:))
    % Nothing fit well! Output a NaN:
    warning('None of the ARCH or GARCH models could be fit.')
    out = NaN;
    return
end

% ------------------------------------------------------------------------------
%% Statistics on retrieved model summaries
% ------------------------------------------------------------------------------
% 'whole things'
out.minLLF = min(LLFs(:));
out.maxLLF = max(LLFs(:));
out.meanLLF = nanmean(LLFs(:));
out.minBIC = min(BICs(:));
out.maxBIC = max(BICs(:));
out.meanBIC = nanmean(BICs(:));
out.minAIC = min(AICs(:));
out.maxAIC = max(AICs(:));
out.meanAIC = nanmean(AICs(:));
out.minK = min(Ks(:));
out.maxK = max(Ks(:));
out.meanK = nanmean(Ks(:));
out.min_meanarchps = min(meanarchps(:));
out.max_meanarchps = max(meanarchps(:));
out.mean_meanarchps = nanmean(meanarchps(:));
out.min_maxarchps = min(maxarchps(:));
out.max_maxarchps = max(maxarchps(:));
out.mean_maxarchps = nanmean(maxarchps(:));
out.min_meanlbqps = min(meanlbqps(:));
out.max_meanlbqps = max(meanlbqps(:));
out.mean_meanlbqps = nanmean(meanlbqps(:));
out.min_maxlbqps = min(maxlbqps(:));
out.max_maxlbqps = max(maxlbqps(:));
out.mean_maxlbqps = nanmean(maxlbqps(:));

% 'bests' (orders)
[a, b] = find(LLFs == min(LLFs(:)),1,'first');
out.bestpLLF = pr(a);
out.bestqLLF = qr(b);

[a, b] = find(AICs == min(AICs(:)),1,'first');
out.bestpAIC = pr(a);
out.bestqAIC = qr(b);

[a, b] = find(BICs == min(BICs(:)),1,'first');
out.bestpAIC = pr(a);
out.bestqAIC = qr(b);

% 'trends' in each direction -- i.e., how much effect on the statistics
% does changing the order (either p or q) have. Sometimes changing q will
% have negligible effect -- we want to quantify this.

out.Ks_vary_p = nanmean(nanstd(Ks)); % mean variation along p direction
out.Ks_vary_q = nanmean(nanstd(Ks')); % mean variation along q direction

end
