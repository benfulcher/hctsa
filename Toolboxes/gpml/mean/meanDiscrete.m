function [m,dm] = meanDiscrete(s, hyp, x)

% Mean function for discrete inputs x. Given a function defined on the
% integers 1,2,3,..,s, the mean function is parametrized as:
%
% m(x) = mu_x,
%
% where mu is a fixed vector of length s.
%
% This implementation assumes that the inputs x are given as integers
% between 1 and s, which simply index the provided vector.
%
% The hyperparameters are:
%
% hyp = [ mu_1
%         mu_2
%         ..
%         mu_s ]
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2016-04-16.
%
% See also covDiscrete.m, meanFunctions.m.

if nargin==0, error('s must be specified.'), end           % check for dimension
if nargin<=2, m = num2str(s); return; end     % report number of hyperparameters
mu = hyp(:); m = mu(x(:));                                       % evaluate mean
dm = @(q) dirder(q,s,x);

function dmdhyp = dirder(q,s,x)
  dmdhyp = zeros(s,1);
  for i=1:s, dmdhyp(i) = sum(q(x==i)); end