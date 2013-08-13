% function ent = MS_shannon(z,bin,depth)
%
% calculate the approximate shannon entropy of a time series using a bin
% bin encoding and depth symbol symbol sequences. The entropy is
% given by 
%    sum P log P
% where P is the probability of any given symbol sequence occuring, 
% and the summation is over all symbol sequences (of depth symbols
% from a bin symbol alphabet).
%
% The binining is a uniform population binning.
%
% defualts:
% bin = 2
% depth = 3
%
% implementation is shannon.c
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 8/10/04
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.