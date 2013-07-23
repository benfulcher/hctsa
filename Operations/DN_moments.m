% DN_Moments
% 
% Output is the moment of the distribution of the input time series.
% Normalizes by the standard deviation
% Uses the moment function in Matlab's Statistics Toolbox
% 
% INPUTS:
% y, the input time series
% n, the moment to calculate (a scalar)
% 

function out = DN_Moments(y,n)
% Ben Fulcher, 2009

out = moment(y,n)/std(y);

end