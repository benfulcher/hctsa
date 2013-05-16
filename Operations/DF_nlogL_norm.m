function nlogL = DF_nlogL_norm(y)
% negative log likelihood of data coming from normal distribution
% Ben Fulcher 2009

[muhat sigmahat] = normfit(y);
nlogL = normlike([muhat sigmahat],y);


end