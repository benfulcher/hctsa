% ST_TrimmedMean
% 
% Outputs the mean of the trimmed time series using the Matlab function trimmean.
% 
% INPUTS:
% y, the input time series
% n, the percent of highest and lowest values in y to exclude from the mean
%     calculation

function out = ST_TrimmedMean(y,n)
% Ben Fulcher, 2008

out = trimmean(y,n);

end