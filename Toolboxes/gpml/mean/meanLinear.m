function [m,dm] = meanLinear(hyp, x)

% Linear mean function. The mean function is parameterized as:
%
% m(x) = sum_i c_i * x_i;
%
% The hyperparameter is:
%
% hyp = [ c_1
%         c_2
%         ..
%         c_D ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<2, m = 'D'; return; end             % report number of hyperparameters
[n,D] = size(x);
if any(size(hyp)~=[D,1]), error('Exactly D hyperparameters needed.'), end
m = x*hyp(:);                                                    % evaluate mean
dm = @(q) x'*q(:);                                      % directional derivative
