function out = DN_quanp(x,p)
% Returns the proportion of the input sequence, x
% lying within p standard deviations of the mean
% Ben Fulcher, 2009

mu = mean(x); % mean of the time series
sig = std(x); % standard deviation of the time series
N = length(x); % length of the time series

out = sum(x > mu-p*sig & x < mu+p*sig)/N;

end