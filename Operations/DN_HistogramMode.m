% DN_HistogramMode
% 
% Measures the mode of the time series using histograms a given numbers
% of bins.
% 
% INPUTS:
% 
% y, the input time series
% 
% nbins, the number of bins to use in the histogram.
% 

function out = DN_HistogramMode(y,nbins)
% Ben Fulcher, October 2009

[dny, dnx] = hist(y,nbins);

% mean of position of maximums (if multiple):
out = mean(dnx(dny == max(dny)));

end