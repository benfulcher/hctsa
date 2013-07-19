% EN_RM_entropy
% 
% Measures the entropy of the time series using a function by Rudy Moddemeijer,
% obtained from http://www.cs.rug.nl/~rudy/matlab/
% This website has code and documentation

function out = EN_RM_entropy(x)
% Wrapper for R. Moddemeijer's entropy code

out = RM_entropy(x);

end