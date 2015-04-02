function [lp,dlp] = priorGauss(mu,s2,x)

% Univariate Gaussian hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(-(x-mu)^2/(2*s2)) / sqrt(2*pi*s2), where
%
% mu(1x1) is the mean parameter, s2(1x1) is the variance parameter and
% x(1xN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-09.
%
% See also PRIORDISTRIBUTIONS.M.

if nargin<2, error('mu and s2 parameters need to be provided'), end
if ~(isscalar(mu)&&isscalar(s2))
  error('mu and s2 parameters need to be scalars')
end
if nargin<3, lp = sqrt(s2)*randn+mu; return, end % return a sample
lp  = -(x-mu).^2/(2*s2) - log(2*pi*s2)/2;
dlp = -(x-mu)/s2;