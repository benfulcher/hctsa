% ST_HighLowMu
% 
% Calculates a statistic related to the mean of the time series data that
% is above the (global) time-series mean compared to the mean of the data that
% is below the global time-series mean.
% 

function out = ST_HighLowMu(y)
% Ben Fulcher, 2008

mu = mean(y);
mhi = mean(y(y > mu));
mlo = mean(y(y < mu));
out = mhi-mu/(mu-mlo);

end