function [m,dm] = meanWarp(mean, g, dg, hyp, x)

% meanWarp - compose a mean function by warping another mean function m0.
%
% m(x) = g( m_0(x) ).
%
% The hyperparameter is:
%
% hyp = [ hyp_m0 ]
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean function to do the actual work.
%
% Copyright (c) by William Herlands & Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<5, m = feval(mean{:}); return, end       % report number of parameters

[m0,dm0] = feval(mean{:},hyp,x);                                   % evaluate m0
m = g(m0);                                                                % mean
dm = @(q) dm0(dg(m0).*q);                               % directional derivative