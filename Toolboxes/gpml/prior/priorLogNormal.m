function [lp,dlp]=priorLogNormal(mu,s2,x)

% Univariate Log-normal hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(-(log(x)-mu)^2/(2*s2)) / (sqrt(2*pi*s2)*x)
%
% where mu(1x1) is the log scale parameter, s2(1x1) is the shape parameter
% and x(1xN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also priorDistributions.m.

if nargin<2, error('mu and s2 parameters need to be provided'), end
if ~(isscalar(mu)&&isscalar(s2))
  error('mu and s2 parameters need to be scalars'), end
if nargin<3, lp = exp(sqrt(s2)*randn+mu); return, end          % return a sample

lx = log(x);
lp  = -(lx-mu).^2/(2*s2) - log(2*pi*s2)/2 - lx;
dlp = -(lx-mu+s2)./(x*s2);
lp(x<0) = -inf; dlp(x<0) = 0;