function A = meanDiscrete(s, hyp, x, i)

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
% Copyright (c) by Roman Garnett, 2014-08-14.
%
% See also COVDISCRETE.M, MEANFUNCTIONS.M.

if nargin==0, error('s must be specified.'), end           % check for dimension
if nargin<=2, A = num2str(s); return; end     % report number of hyperparameters
mu = hyp(:);
if nargin==3
  A = mu(x(:));                                                  % evaluate mean
else
  A = zeros(numel(x),1);                                            % derivative
  A(x==i) = 1;
end