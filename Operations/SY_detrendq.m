function out = SY_detrendq(x)
% Returns the standard deviation of the linearly-detrended time series
% as a function of its original standard deviation
% Ben Fulcher, 2009

out = std(detrend(x)) / std(x);
    
end
