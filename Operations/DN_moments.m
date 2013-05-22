function out = DN_moments(y,n)
% Returns moments of the distribution of values in the input time series, y
% Uses the moment function in Matlab's Statistics Toolbox
% Normalizes by the standard deviation
% Ben Fulcher, 2009

out = moment(y,n)/std(y);

end