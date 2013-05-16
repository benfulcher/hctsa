function out = SY_detrendq(x)
    out = std(detrend(x))/std(x);
end
