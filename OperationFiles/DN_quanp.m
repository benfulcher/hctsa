function out = DN_quanp(y,p)

mu = mean(y);
sig = std(y);
% lims=quantile(y,[0.5-th 0.5+th]);
a = find(y>mu-p*sig & y<mu+p*sig);
out = length(a)/length(y); % proportion of y lying within p stds of the mean

end