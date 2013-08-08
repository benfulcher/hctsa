% EN_RM_entropy
% 
% Measures the entropy of the time series using a function by Rudy Moddemeijer
% 
% Original code, now RM_entropy, was obtained from:
% http://www.cs.rug.nl/~rudy/matlab/
% 
% The above website has code and documentation for the function.
% 
% INPUTS:
% y, the input time series
% 

function out = EN_RM_entropy(y)
% Wrapper for R. Moddemeijer's entropy code, RM_entropy

out = RM_entropy(y);

end