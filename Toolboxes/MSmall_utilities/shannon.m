% function ent=shannon(z,bin,depth)
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
% 8/10/04
