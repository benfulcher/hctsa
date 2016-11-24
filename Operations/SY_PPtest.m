function out = SY_PPtest(y,lags,model,testStatistic)
% SY_PPtest   Phillips-Peron unit root test.
%
% Uses the pptest code from Matlab's Econometrics Toolbox.
%
%---INPUTS:
% y, the input time series
%
% lags, a vector of lags
%
% model, a specified model:
%               'ar': autoregressive
%               'ard': autoregressive with drift, or
%               'ts': trend stationary,
%               (see Matlab documentation for information)
%
% testStatistic, the test statistic:
%               't1': the standard t-statistic, or
%               't2' a lag-adjusted, 'unStudentized' t statistic.
%               (see Matlab documentation for information)
%
%---OUTPUTS: statistics on the p-values and lags obtained from the set of tests, as
% well as measures of the regression statistics.

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
warning('off','econ:pptest:LeftTailStatTooSmall')
[h, pValue, stat, ~, reg] = pptest(y,'lags',lags,'model',model,'test',testStatistic);
warning('on','econ:pptest:LeftTailStatTooSmall')

% ------------------------------------------------------------------------------
%% Get outputs
% ------------------------------------------------------------------------------
nout = length(h);

if nout == 1
    % Just return the results from this single test
    out.pvalue = pValue;
    out.stat = stat;
    out.coeff1 = reg.coeff(1); % could be multiple, depending on the model
    out.loglikelihood = reg.LL;
    out.AIC = reg.AIC;
    out.BIC = reg.BIC;
    out.HQC = reg.HQC;
    out.rmse = reg.RMSE;

else
    % Return statistics on the set of outputs
    out.maxpValue = max(pValue);
    out.minpValue = min(pValue);
    out.meanpValue = mean(pValue);
    out.stdpValue = std(pValue);
    imaxp = find(pValue == max(pValue),1,'first');
    iminp = find(pValue == min(pValue),1,'first');
    out.lagmaxp = lags(imaxp);
    out.lagminp = lags(iminp);

    out.meanstat = mean(stat);
    out.maxstat = max(stat);
    out.minstat = min(stat);

    % Some regression statistics
    % These are all highly correlated; hctsa library records just minBIC by default
    out.meanloglikelihood = mean(vertcat(reg.LL));
    out.minAIC = min(vertcat(reg.AIC));
    out.minBIC = min(vertcat(reg.BIC));
    out.minHQC = min(vertcat(reg.HQC));

    % Somehow these are highly correlated, only minrmse is included in hctsa
    % library by default
    out.minrmse = min(vertcat(reg.RMSE));
    out.maxrmse = max(vertcat(reg.RMSE));
end

end
