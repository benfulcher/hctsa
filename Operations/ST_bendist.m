function out = ST_bendist(y)
% Returns a statistics related to the mean of the time series above 
% and below its mean
% Ben Fulcher, 2008

mu = mean(y);
mhi = mean(y(y > mu));
mlo = mean(y(y < mu));
out = mhi-mu/(mu-mlo);

end