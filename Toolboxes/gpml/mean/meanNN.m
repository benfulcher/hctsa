function [m,dm] = meanNN(c,v, hyp, x)

% Nearest neighbor mean function. The mean function is parameterized as:
%
% m(z) = v_j, j = arg min_i d(ci,x) where d is the Euclidean distance and ci is
%                                   the ith cluster center.
%
% The hyperparameters are:
%
% hyp = [ ]
%
% Copyright (c) by Hannes Nickisch, 2016-04-16.
%
% See also meanFunctions.m.

if nargin<4, m = '0'; return; end             % report number of hyperparameters 
if numel(hyp)~=0, error('No hyperparameters needed for this model.'), end

[junk,j] = min(sq_dist(c',x')); m = v(j); m = m(:);              % evaluate mean
dm = @(q) zeros(0,1);                                   % directional derivative
