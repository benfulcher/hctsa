function out = ST_mode(y,nbins)
% Finds the mode of the binned distribution with nbins specified
% Updated Ben Fulcher October 2009

[dny, dnx] = hist(y,nbins);
out = mean(dnx(dny == max(dny)));

end