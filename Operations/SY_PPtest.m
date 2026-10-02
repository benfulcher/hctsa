function out = SY_PPtest(y, lags, model, testStatistic)
% SY_PPtest   Phillips-Perron unit root test.
%
% Uses the pptest function from MATLAB's Econometrics Toolbox. The null hypothesis
% is that the series is a unit-root process; the alternative is that it is
% stationary (small p-values reject a unit root). The test is run for each of a
% set of numbers of autocovariance lags, and the outputs summarize the set of
% tests. (With a single lag, the results of that one test are returned instead.)
%
% ---INPUTS:
% y, the input time series
%
% lags, a vector of lags: the numbers of autocovariance lags included in the
%       Newey-West estimator of the long-run variance (default: 0:5)
%
% model, a specified model (default: 'ar'):
%               'ar': autoregressive
%               'ard': autoregressive with drift, or
%               'ts': trend stationary,
%               (see MATLAB documentation for information)
%
% testStatistic, the test statistic (default: 't1'):
%               't1': the standard t-statistic, or
%               't2': a lag-adjusted, 'unStudentized' t statistic.
%               (see MATLAB documentation for information)
%
% ---OUTPUTS:
% For a vector of lags:
% minpValue, meanpValue: the minimum and mean p-value across the tests
% lagmaxp, lagminp: the lag at which the p-value is largest and smallest (the
%       first such lag in the case of a tie)
% meanstat: the mean test statistic across the tests
% minBIC: the minimum, across the tests, of the Bayesian information criterion of
%       the test regression, per observation
% For a single lag:
% pvalue, stat: the p-value and test statistic
% coeff1: the first regression coefficient
% loglikelihood, AIC, BIC, HQC: the log likelihood and information criteria of the
%       test regression, per observation
% rmse: the root-mean-square error of the regression

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
%% Check that an Econometrics Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('econometrics_toolbox');

% ------------------------------------------------------------------------------
%% Inputs
% ------------------------------------------------------------------------------
% The number of autocovariance lags to include in the Newey-West estimator
% of the long-run variance, lags.
if nargin < 2 || isempty(lags)
	lags = (0:5); % use 5 autoregressive lags.
end

% The model variant, model. Can be 'ar': autoregressive, 'ard':
% autoregressive with drift, or 'ts': trend stationary.
if nargin < 3 || isempty(model)
	model = 'ar'; % autoregressive
end

% The test statistics, testStatistic. Can be 't1': standard t statistics; or
% 't2': a lag-adjusted, 'unstudentized' t statistic.
if nargin < 4 || isempty(testStatistic)
	testStatistic = 't1'; % standard t statistic
end

% ------------------------------------------------------------------------------
%% Run the test
% ------------------------------------------------------------------------------
warning('off', 'econ:pptest:LeftTailStatTooSmall')
[h, pValue, stat, ~, reg] = pptest(y, 'lags', lags, 'model', model, 'test', testStatistic);
warning('on', 'econ:pptest:LeftTailStatTooSmall')

% ------------------------------------------------------------------------------
%% Get outputs
% ------------------------------------------------------------------------------
nout = length(h);

if nout == 1
	% Just return the results from this single test
	out.pvalue = pValue;
	out.stat = stat;
	out.coeff1 = reg.coeff(1); % could be multiple, depending on the model
	% Log-likelihood and the information criteria are extensive: they are sums
	% over observations, so they grow in direct proportion to the time-series
	% length regardless of how well the model fits (minBIC measured eta^2 = 0.973
	% against N across series of different lengths). Reported per observation, which
	% is the standard intensive form and the quantity model comparison actually
	% depends on.
	numObs = length(y);
	out.loglikelihood = reg.LL / numObs;
	out.AIC = reg.AIC / numObs;
	out.BIC = reg.BIC / numObs;
	out.HQC = reg.HQC / numObs;
	out.rmse = reg.RMSE;

else
	% Return statistics on the set of outputs. maxpValue/stdpValue dropped
	% (r >= 0.98 with meanpValue on real-world series); maxstat/minstat dropped
	% (r >= 0.97 with meanstat).
	out.minpValue = min(pValue);
	out.meanpValue = mean(pValue);
	imaxp = find(pValue == max(pValue), 1, 'first');
	iminp = find(pValue == min(pValue), 1, 'first');
	out.lagmaxp = lags(imaxp);
	out.lagminp = lags(iminp);

	out.meanstat = mean(stat);

	% Regression statistics: meanloglikelihood/minAIC/minHQC/minrmse/maxrmse
	% dropped -- confirmed r >= 0.998 with minBIC (and with each other) on
	% real-world series, matching this function's own longstanding comment that
	% these are all highly correlated. Per observation -- see the note in
	% the single-test branch above.
	numObs = length(y);
	out.minBIC = min(vertcat(reg.BIC)) / numObs;
end

end
