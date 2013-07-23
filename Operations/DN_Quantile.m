% DN_Quantile
% 
% Calculates the quantile value at a specified proportion p using
% the Statistics Toolbox function, quantile.
% 
% INPUTS:
% y, the input time series
% p, the quantile proportion
% 

function out = DN_Quantile(y,p)
% Ben Fulcher 5/8/09

if nargin < 2
    fprintf(1,'Using quantile p = 0.5 (median) by default\n');
    p = 0.5;
end
if (p < 0) || (p > 1)
    error('p must specify a proportion, 0 <= p <= 1');
end

out = quantile(y,p);

end