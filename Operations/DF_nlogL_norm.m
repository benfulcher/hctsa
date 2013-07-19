% DF_nlogL_norm
% 
% Fits a Gaussian distribution to the data using the normfit function in
% MATLAB's Statistics Toolbox and returns the negative log likelihood of the
% data coming from a Gaussian distribution using the normlike function.
% 
% INPUT: y, the time series.

function nlogL = DF_nlogL_norm(y)
% Ben Fulcher, 2009

[muhat, sigmahat] = normfit(y);
nlogL = normlike([muhat, sigmahat],y);

end