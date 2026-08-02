function [lp,dlp] = priorWeibull(lam,k,x)

% Univariate Weibull hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = (k/lam) * (x/lam)^(k-1) * exp(-(x/lam)^k)
%
% where lam(1x1) is the scale parameter, k(1x1) is the scale parameter
% and x(1xN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also priorDistributions.m.

if nargin<2, error('lam and k parameters need to be provided'), end
if ~(isscalar(lam)&&isscalar(k))
  error('lam and k parameters need to be scalars'), end
if nargin<3,lp = lam*(-log(rand))^(1/k); return, end           % return a sample

lp = log(k/lam) + (k-1)*log(x/lam) - (x/lam).^k;
dlp = (k-1)./x - (k/lam)*(x/lam).^(k-1);
lp(x<0) = -inf; dlp(x<0) = 0;