function [m,dm] = meanZero(hyp, x)
%MEANZERO Zero mean function.
%
% Report number of hyperparameters
%  s = MEANZERO ()
%  s = MEANZERO (hyp)
%
% Evaluation of m(x)
%  m = MEANZERO (hyp, x)
%
% Evaluation of function derivatives dm (w.r.t. hyp)
%    [m, dm] = MEANZERO (hyp, x)
%
% Call meanFunctions to get an explanation of outputs in each mode.
%
% This mean function does not have any parameters.
%
% m(x) = 0
%
% See also MEANFUNCTIONS

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-15.

if nargin<2, m = '0'; return; end             % report number of hyperparameters
m = zeros(size(x,1),1);                                                   % mean
dm = @(q) zeros(0,1);                                   % directional derivative
