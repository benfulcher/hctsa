% DN_pleft
% 
% Measures the maximum distance from the mean at which a given fixed proportion,
% p, of the time-series data points are further.
% Normalizes by the standard deviation of the time series
% (could generalize to separate positive and negative deviations in future)
% Uses the quantile function from Matlab's Statistics Toolbox
%
% INPUTS:
% y, the input time series
% th, the proportion of data further than p from the mean
%           (output p, normalized by standard deviation)

function out = DN_pleft(y,th)
% Ben Fulcher, 2009

if nargin < 2 || isempty(th)
    th = 0.1; % default
end

p = quantile(abs(y-mean(y)),1-th);

% A proportion, th, of the data lie further than p from the mean
out = p/std(y);

end