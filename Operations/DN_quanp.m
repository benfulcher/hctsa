function out = DN_quanp(y,p)
% Returns the proportion of y lying within p standard deviations of the mean
% Ben Fulcher, 2009

mu = mean(y);
sig = std(y);
N = length(y);

out = sum(y>mu-p*sig & y<mu+p*sig)/N;


% lims = quantile(y,[0.5-th 0.5+th]);
% a = find(y>mu-p*sig & y<mu+p*sig);

end