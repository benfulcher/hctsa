function out = ST_bendist(y)

mu = mean(y);
mhi = mean(y(y>mu));
mlo = mean(y(y<mu));
out = mhi-mu/(mu-mlo);

end