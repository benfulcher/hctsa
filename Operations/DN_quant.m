function out = DN_quant(y,p)
% Calculates the quantile, pure and simple
% Uses the quantile function from Matlab's Statistics Toolbox
% Ben Fulcher 5/8/09

out = quantile(y,p);

end