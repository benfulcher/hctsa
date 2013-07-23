% EN_PermEn
% 
% Computes the permutation entropy of order, ord, of a time series.
% 
% "Permutation Entropy: A Natural Complexity Measure for Time Series"
% C. Bandt and B. Pompe, Phys. Rev. Lett. 88(17) 174102 (2002)
% 
% Code is adapted from logisticPE.m code obtained from
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m
% 
% INPUTS:
% y, a time series
% ord, the order of permutation entropy
% 

function permen = EN_PermEn(y,ord)
% Ben Fulcher, 2009

% Ensure y is a column vector
if size(y,1) > size(y,2);
    y = y';
end

% Use the Bruce Land and Damian Elias code to calculate the permutation entropy:
permen = LA_permen(y);

end