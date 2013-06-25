function nlogL = DF_nlogL_norm(y)
% Calculates the negative log likelihood of data coming from normal distribution
% Uses the normfit and normlike functions from Matlab's Statistics Toolbox
% Ben Fulcher 2009

[muhat, sigmahat] = normfit(y);
nlogL = normlike([muhat, sigmahat],y);

end