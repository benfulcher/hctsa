function out = ST_cv(x,k)
% Coefficient of Variation, CV
% k a positive integer: the 'moment' of cv

out = (std(x))^k / (mean(x))^k;

end