% DN_quant
% 
% Calculates the quantile value at a specified proportion p using
% the Statistics Toolbox function, quantile.

function out = DN_quant(y,p)
% Ben Fulcher 5/8/09

if nargin < 2
    fprintf(1,'Using quantile p = 0.5 (median) by default\n');
    p = 0.5;
end

out = quantile(y,p);

end