function out = DN_quant(y,p)
% Outputs the quantile of the set of values in y
% Uses the quantile function from Matlab's Statistics Toolbox
% Ben Fulcher 5/8/09

out = quantile(y,p);

end