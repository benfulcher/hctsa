% cmp = MS_complexitybs(x,n);
%
% calculate the Lempel-Ziv complexity of the n-symbol stream x. 
%
% x\in{0,1,2,...,(n-1)}
%
% cmp is the normalised complexity, that is the number of distinct
% symbol sequences in x, divided by the expected number of distinct 
% symbols for a noise sequence.
%
% Algorithm is implemented in MS_complexitybs.c
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 24/9/03
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.