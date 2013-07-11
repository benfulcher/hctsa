function out = DN_quanp(x,p)
% Returns the proportion of the input sequence, x
% lying within p standard deviations of the mean
% Ben Fulcher, 2009

mu = mean(x);
sig = std(x);
N = length(x);

out = sum(x > mu-p*sig & x < mu+p*sig)/N;

end