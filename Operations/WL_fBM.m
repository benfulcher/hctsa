% WL_fBM
% 
% Uses the wfbmesti function from Matlab's Wavelet Toolbox to estimate the
% parameters of fractional Gaussian Noise, or fractional Brownian motion in a
% time series.
% 
% INPUT:
% y, the time series to analyze.
% 
% All three outputs of wfbmesti are returned from this function.
% 

function out = WL_fBM(y)
% Relies on the Wavelet Toolbox function, wfbmesti
% Ben Fulcher 23/1/2010

% Parameter estimation of fractional Brownian motion
hest = wfbmesti(y);
out.p1 = hest(1);
out.p2 = hest(2);
out.p3 = hest(3);

end