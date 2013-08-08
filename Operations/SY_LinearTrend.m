% SY_LinearTrend
% 
% Linearly detrends the time series using the Matlab algorithm detrend,
% and returns the ratio of standard deviations before and after the linear
% detrending.
% 
% If a strong linear trend is present in the time series, this  operation should
% output a low value.
% 
% INPUT:
% x, the input time series
% 

function out = SY_LinearTrend(x)
% Ben Fulcher, 2009

out = std(detrend(x)) / std(x);
    
end
