function out = EC_kpsstest(y,lags)
% Uses the MATLAB Econometrics Toolbox to perform the kpss stationarity
% test. Null hypothesis is that that univariate time series given by y is
% trend stationary; against the alternative that it's a nonstationary,
% unit-root process.
% Ben Fulcher 26/2/2010

% lags, can be a scalar, which returns the values of that test as output;
% or a vector, which returns statistics on the outputs as applied across
% this range of lags.
if nargin<2 || isempty(lags)
    lags = 0;
end

%% (1) Perform the test(s)
[h,pValue,stat,cValue] = kpsstest(y,'lags',lags);


%% (2) Return statistics on outputs of test(s)
if length(lags)>1
    % return statistics on outputs
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

% keyboard


end