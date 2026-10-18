function out = MF_GARCHcompare(y, preProc, pr, qr, randomSeed, beVocal)
% MF_GARCHcompare   How well GARCH models of different orders describe the changing variance of the series.
%
% Fits a set of GARCH(p,q) models, with zero mean and Gaussian innovations, to the
% time series, and returns statistics on the goodness of fit across a range of p (the
% number of lagged variances) and q (the number of lagged squared innovations). All
% combinations of p and q are compared, and the outputs are summaries across the
% grid of fitted models, and the orders that fit best.
%
% The models are fitted by maximum likelihood with BF_GARCHFit (see MF_GARCHfit), in order
% of increasing p and q, each also started from the fits of its smaller neighbors (nested
% models), so the log-likelihood never decreases as an order increases. The
% GARCH(min(p),min(q)) model is tested against constant variance as in MF_GARCHfit; if it does
% not beat it (hasStructure is 0), the constant term and the residual statistics are those of
% constant variance, since the coefficients of the models are not identified.
%
% This operation focuses on the GARCH/variance component. It therefore optionally
% pre-whitens the time series (BF_Whiten, which detrends it first) and z-scores it,
% so that a constant-mean process is assumed.
%
% ---INPUTS:
% y, the input time series
%
% preProc, a preprocessing to apply:
%           (i) 'none': no preprocessing is performed, and the series is not
%               detrended (default),
%           (ii) 'ar': detrends the series, then performs a preprocessing that
%               maximizes AR(2) whiteness.
%
% pr, a vector of model orders, p, to compare (default 1:3)
%
% qr, a vector of model orders, q, to compare (default 1:3)
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed (used
%           in the whitening, BF_Whiten)
%
% beVocal, whether to print messages about failed fits to the command line
%           (default 0)
%
% ---OUTPUTS: statistics across the (p,q) models that fitted successfully. The
% log-likelihood, AIC and BIC are per observation, i.e., divided by the series length.
% hasStructure: 1 if the GARCH(min(p),min(q)) model beats constant variance, otherwise 0
% minLLF, maxLLF, meanLLF: the minimum, maximum and mean log-likelihood
% minBIC, maxBIC, meanBIC: the minimum, maximum and mean Bayesian information
%       criterion (BIC)
% minAIC, maxAIC, meanAIC: the minimum, maximum and mean Akaike information
%       criterion (AIC)
% minK, maxK, meanK: the minimum, maximum and mean constant term of the variance
%       equation
% min_meanarchps, max_meanarchps, mean_meanarchps: the minimum, maximum and mean, across
%       models, of the mean p-value (over lags 1 to 20) of Engle's ARCH test on the
%       standardized residuals
% min_maxarchps, max_maxarchps, mean_maxarchps: the same for the maximum p-value over
%       the 20 lags
% min_meanlbqps, max_meanlbqps, mean_meanlbqps: the minimum, maximum and mean, across
%       models, of the mean p-value (over lags 1 to 20) of the Ljung-Box Q-test on the
%       squared standardized residuals
% min_maxlbqps, max_maxlbqps, mean_maxlbqps: the same for the maximum p-value over the
%       20 lags
% bestpAIC, bestqAIC: the orders p and q of the model with the minimum AIC
% bestpBIC, bestqBIC: the orders p and q of the model with the minimum BIC
% Ks_vary_p, Ks_vary_q: how much the constant term varies with p, and with q: the
%       standard deviation of the constant across one order, averaged over the other
%
% ---NOTES:
% The orders chosen by AIC and BIC are among the grid and the constant-variance model
% (orders 0 and 0), which is chosen when the GARCH models do not improve on it enough to
% pay for their parameters; ties go to the smaller model.
%
% Registered in hctsa's feature set: maxLLF, bestpAIC, bestqAIC, bestpBIC, bestqBIC,
% min_meanarchps and max_meanarchps. The constant terms (minK, maxK, meanK) follow from the
% persistence, and the other summaries repeat these (FeatureSets/renamed_operations.csv gives
% the reason for each).

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

pr = sort(pr(:)'); % fitted in order of increasing p and q
qr = sort(qr(:)');

% ------------------------------------------------------------------------------
%% (1) Data preprocessing
% ------------------------------------------------------------------------------
y = BF_Whiten(y, preProc, 0, randomSeed); % Use a basic preprocessing to whiten the time series
y = zscore(y); % make sure the (whitened) time series is z-scored
N = length(y); % could be different to original (e.g., if chose a differencing above)

% Now have the preprocessed time series saved over y.

% ------------------------------------------------------------------------------
%% (2) Fit the GARCH(p,q) models
% ------------------------------------------------------------------------------
D = BF_GARCHFit('data', y);
nullFit = BF_GARCHFit('null', D, false); % constant variance

np = length(pr);
nq = length(qr);
fits = cell(np, nq);
for i = 1:np
	p = pr(i); % garch order
	for j = 1:nq
		q = qr(j); % arch order
		try
			% Start from the fits of the models one order smaller (nested in this one)
			S = BF_GARCHFit('spec', D, p, q, false, false);
			warm = zeros(S.k, 0);
			if j > 1 && ~isempty(fits{i, j-1})
				warm(:, end+1) = BF_GARCHFit('embed', fits{i, j-1}.theta, fits{i, j-1}.spec, S);
			end
			if i > 1 && ~isempty(fits{i-1, j})
				warm(:, end+1) = BF_GARCHFit('embed', fits{i-1, j}.theta, fits{i-1, j}.spec, S);
			end
			fits{i, j} = BF_GARCHFit('mle', D, p, q, false, false, warm, (i == 1 && j == 1));
		catch
			if beVocal
				fprintf(1, 'Bad fit at p = %u, q = %u\n', p, q);
			end
		end
	end
end

% Check if fits went well...
if isempty(fits{1, 1})
	% Nothing fit well! Output a NaN:
	warning('None of the ARCH or GARCH models could be fit.')
	out = NaN;
	return
end

% Is there any conditional heteroskedasticity to model? If not, the constant and the
% residual statistics are those of constant variance
[~, hasStructure] = BF_GARCHFit('structure', fits{1, 1}, nullFit, D, false);

% ------------------------------------------------------------------------------
%% (3) Store derived statistics on the fitted models
% ------------------------------------------------------------------------------
LLFs = NaN(np, nq); % log-likelihood
AICs = NaN(np, nq); % AIC
BICs = NaN(np, nq); % BIC
Ks = NaN(np, nq); % constant term in variance
meanarchps = NaN(np, nq); % mean over 20 lags from Engle's ARCH test on standardized innovations
maxarchps = NaN(np, nq); % maximum p-value over 20 lags from Engle's ARCH test on standardized innovations
meanlbqps = NaN(np, nq); % mean lbq p-value over 20 lags from Q-test on squared standardized innovations
maxlbqps = NaN(np, nq); % maximum p-value over 20 lags from Q-test on squared standardized innovations

for i = 1:np
	for j = 1:nq
		ml = fits{i, j};
		if isempty(ml)
			continue % didn't fit successfully; everything stays NaN
		end
		% (the likelihood and information criteria are those of the maximum-likelihood
		% fit; the constant and the residual tests those of the reported fit)
		if hasStructure
			fit = ml;
		else
			fit = nullFit;
		end
		LLFs(i, j) = ml.llf;
		AICs(i, j) = -2*ml.llf + 2*ml.k;
		BICs(i, j) = -2*ml.llf + ml.k*log(N);
		Ks(i, j) = fit.theta(1);

		% Summaries of the standardized errors (the absolute values are less
		% meaningful than the comparison between models)
		stde = y ./ sqrt(fit.h); % standardized residuals
		pArch = BF_ARCHTest(stde, 20); % Engle's ARCH test, lags 1:20
		pLbq = BF_LjungBoxTest(stde.^2, 20); % Ljung-Box Q-test of the squared residuals, lags 1:20
		meanarchps(i, j) = mean(pArch);
		maxarchps(i, j) = max(pArch);
		meanlbqps(i, j) = mean(pLbq);
		maxlbqps(i, j) = max(pLbq);
	end
end

% ------------------------------------------------------------------------------
%% Statistics on retrieved model summaries
% ------------------------------------------------------------------------------
% Log-likelihoods and information criteria are sums over observations, so they
% scale with the series length; express them per observation so they are
% comparable across lengths (order selection below is unaffected)
LLFn = LLFs / N;
AICn = AICs / N;
BICn = BICs / N;

out.hasStructure = double(hasStructure);

% 'whole things'
out.minLLF = min(LLFn(:));
out.maxLLF = max(LLFn(:));
out.meanLLF = mean(LLFn(:),'omitnan');
out.minBIC = min(BICn(:));
out.maxBIC = max(BICn(:));
out.meanBIC = mean(BICn(:),'omitnan');
out.minAIC = min(AICn(:));
out.maxAIC = max(AICn(:));
out.meanAIC = mean(AICn(:),'omitnan');
out.minK = min(Ks(:));
out.maxK = max(Ks(:));
out.meanK = mean(Ks(:),'omitnan');
out.min_meanarchps = min(meanarchps(:));
out.max_meanarchps = max(meanarchps(:));
out.mean_meanarchps = mean(meanarchps(:),'omitnan');
out.min_maxarchps = min(maxarchps(:));
out.max_maxarchps = max(maxarchps(:));
out.mean_maxarchps = mean(maxarchps(:),'omitnan');
out.min_meanlbqps = min(meanlbqps(:));
out.max_meanlbqps = max(meanlbqps(:));
out.mean_meanlbqps = mean(meanlbqps(:),'omitnan');
out.min_maxlbqps = min(maxlbqps(:));
out.max_maxlbqps = max(maxlbqps(:));
out.mean_maxlbqps = mean(maxlbqps(:),'omitnan');

% 'bests' (orders): among the grid and the constant-variance model (orders 0, 0), which
% is chosen when the models do not improve on it enough to pay for their parameters.
% Ties go to the smaller model.
nullAIC = -2*nullFit.llf + 2*nullFit.k;
nullBIC = -2*nullFit.llf + nullFit.k*log(N);
critNames = {'AIC', 'BIC'};
critValues = {AICs, BICs};
critNull = [nullAIC, nullBIC];
for c = 1:2
	best = critNull(c); bp = 0; bq = 0;
	for j = 1:nq
		for i = 1:np
			if critValues{c}(i, j) < best
				best = critValues{c}(i, j); bp = pr(i); bq = qr(j);
			end
		end
	end
	out.(['bestp' critNames{c}]) = bp;
	out.(['bestq' critNames{c}]) = bq;
end

% 'trends' in each direction -- i.e., how much effect on the statistics
% does changing the order (either p or q) have. Sometimes changing q will
% have negligible effect -- we want to quantify this.

out.Ks_vary_p = mean(std(Ks,0,'omitnan'),'omitnan'); % mean variation along p direction
% mean variation along q direction:
out.Ks_vary_q = mean(std(Ks',0,'omitnan'),'omitnan'); %#ok<UDIM>

end
