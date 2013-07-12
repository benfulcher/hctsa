function out = DN_pleft(y,th)
% Returns the distance from the mean at which a proportion p of the time series values lie further from the mean than this point
% Normalizes by the standard deviation of the time series
% (could generalize to separate positive and negative deviations)
% Uses the quantile function from Matlab's Statistics Toolbox
% Ben Fulcher, 2009

mu = mean(y);
sig = std(y);

p = quantile(abs(y-mu),1-th);

% A proportion, th, of the data lie further than p from the mean
out = p/sig;

end