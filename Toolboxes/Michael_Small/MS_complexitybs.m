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
% M. Small
% ensmall@polyu.edu.hk
% 24/9/03
