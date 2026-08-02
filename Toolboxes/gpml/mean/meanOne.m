function [m,dm] = meanOne(hyp, x)

% One mean function. The mean function does not have any parameters.
%
% m(x) = 1
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<2, m = '0'; return; end             % report number of hyperparameters
m = ones(size(x,1),1);                                                    % mean
dm = @(q) zeros(0,1);                                   % directional derivative
