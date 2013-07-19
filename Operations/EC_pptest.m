% EC_pptest
% 
% Performs the Phillips-Peron unit root test for a time series via the code
% pptest from MATLAB's Econometrics Toolbox.
% 
% INPUTS:
% y, the input time series
% lags, a vector of lags
% model, a specified model:
%               'ar': autoregressive
%               'ard': autoregressive with drift, or
%               'ts': trend stationary,
%               (see Matlab documentation for information)
% teststat, the test statistic:
%               't1': the standard t-statistic, or
%               't2' a lag-adjusted, 'unStudentized' t statistic.
%               (see Matlab documentation for information)
%               
% Outputs are statistics on the p-values and lags obtained from the set of tests, as
% well as measures of the regression statistics.

function out = EC_pptest(y,lags,model,teststat)
% Ben Fulcher, 1/3/2010

%% Inputs
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

% The test statistics, teststat. Can be 't1': standard t statistics; or
% 't2': a lag-adjusted, 'unstudentized' t statistic.
if nargin < 4 || isempty(teststat)
    teststat = 't1'; % standard t statistic
end


%% Run the test
[h, pValue, stat, cValue, reg] = pptest(y,'lags',lags,'model',model,'test',teststat);

%% Get outputs
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
    
    % some regression statistics
%     out.meancoeff1 = mean(reg.coeff(1)); % could be multiple, depending on the model
    out.meanloglikelihood = mean(vertcat(reg.LL));
    out.minAIC = min(vertcat(reg.AIC));
    out.minBIC = min(vertcat(reg.BIC));
    out.minHQC = min(vertcat(reg.HQC));
    out.minrmse = min(vertcat(reg.RMSE));
    out.maxrmse = max(vertcat(reg.RMSE));
end

end