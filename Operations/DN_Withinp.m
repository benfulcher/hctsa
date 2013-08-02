% DN_Withinp
% 
% Measures the proportion of the time-series data points that lie within
% p standard deviations of its mean.
% 
% INPUTS:
% x, the input time series
% p, the number (proportion) of standard deviations.
% 

function out = DN_Withinp(x,p)
% Ben Fulcher, 2009

if nargin < 2 || isempty(p)
    p = 1; % 1 std from mean
end

mu = mean(x); % mean of the time series
sig = std(x); % standard deviation of the time series
N = length(x); % length of the time series

out = sum(x >= mu-p*sig & x <= mu+p*sig)/N;

end