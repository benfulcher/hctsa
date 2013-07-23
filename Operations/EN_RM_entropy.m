% EN_RM_entropy
% 
% Measures the entropy of the time series using a function by Rudy Moddemeijer
% 
% Original code, now RM_entropy, was obtained from http://www.cs.rug.nl/~rudy/matlab/
% This website has code and documentation for the function.
% 

function out = EN_RM_entropy(x)
% Wrapper for R. Moddemeijer's entropy code

out = RM_entropy(x);

end