function [m,dm] = meanConst(hyp, x)

% Constant mean function. The mean function is parameterized as:
%
% m(x) = c
%
% The hyperparameter is:
%
% hyp = [ c ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<2, m = '1'; return; end             % report number of hyperparameters 
if numel(hyp)~=1, error('Exactly one hyperparameter needed.'), end
c = hyp;
m = c*ones(size(x,1),1);                                                  % mean
dm = @(q) sum(q);                                       % directional derivative
