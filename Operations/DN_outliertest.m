% DN_OutlierTest
% 
% Removes the p% of highest and lowest values in the time series (i.e., 2*p%
% removed from the time series in total) and returns the ratio of either the
% mean or the standard deviation of the time series, before and after this
% transformation.
% 
% INPUTS:
% y, the input time series
% p, the percentage of values to remove beyond upper and lower percentiles
% 

function out = DN_OutlierTest(y,p)
% Ben Fulcher, 2009

if nargin < 2
    p = 2; % by default, remove 2% of values from upper and lower percentiles
end

% mean of the middle (100-2*p)% of the data
out.mean = mean(y(y > prctile(y,p) & y < prctile(y,100-p)));

% std of the middle (100-2*p)% ofthe data
out.std = std(y(y > prctile(y,p) & y < prctile(y,100-p))) / std(y); % [although std(y) should be 1]

end