% EC_kpsstest
% 
% Performs the KPSS stationarity test, of Kwiatkowski, Phillips, Schmidt, and Shin,
% "Testing the null hypothesis of stationarity against the alternative of a
% unit root: How sure are we that economic time series have a unit root?"
% Kwiatkowski, Denis and Phillips, Peter C. B. and Schmidt, Peter and Shin, Yongcheol
% J. Econometrics, 54(1-3) 159 (2002)
% 
% Uses the function kpsstest from MATLAB's Econometrics Toolbox. The null
% hypothesis is that a univariate time series is trend stationary, the
% alternative hypothesis is that it is a non-stationary unit-root process.
% 
% The code can implemented for a specific time lag, tau. Alternatively, measures
% of change in p-values and test statistics will be outputted if the input is a
% vector of time lags.
% 
% INPUTS:
% y, the input time series
% lags, can be either a scalar (returns basic test statistic and p-value), or
%                   vector (returns statistics on changes across these time lags)

function out = EC_kpsstest(y,lags)
% Ben Fulcher, 26/2/2010

% Check inputs
if nargin < 2 || isempty(lags)
    lags = 0;
end

%% (1) Perform the test(s)
[h, pValue, stat, cValue] = kpsstest(y,'lags',lags);


%% (2) Return statistics on outputs of test(s)
if length(lags) > 1
    % Return statistics on outputs
    out.maxpValue = max(pValue);
    out.minpValue = min(pValue);
    out.maxstat = max(stat);
    out.minstat = min(stat);
    out.lagmaxstat = lags(stat == max(stat)); % lag at max test statistic
    out.lagminstat = lags(stat == min(stat));
else
    % return the statistic and pvalue
    out.stat = stat;
    out.pValue = pValue;
end

end