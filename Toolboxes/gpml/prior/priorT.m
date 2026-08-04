function [lp,dlp] = priorT(mu,s2,nu,x)

% Univariate Student's t hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = ( 1 + (x-mu)^2/(s2*(nu-2)) )^(-(nu+1)/2) / Z, where
%      Z = gamma(nu/2) * sqrt(pi*(nu-2)*s2) / gamma((nu+1)/2),
%
% mu(1x1) is the mean parameter, s2(1x1) is the variance parameter,
% nu(1x1) > 2 is the degrees of freedom parameter and x(1xN) contains query
% hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-10-16.
%
% See also priorDistributions.m.

if nargin<3, error('mu, s2 and nu parameters need to be provided'), end
if ~(isscalar(mu)&&isscalar(s2)&&isscalar(nu))
  error('mu, s2 and nu parameters need to be scalars')
end
if nargin<4, lp = mu + sqrt(s2*(nu-2)/priorGamma(nu/2,2))*randn; return, end

lZ = gammaln((nu+1)/2)-gammaln(nu/2)-log(pi*(nu-2)*s2)/2;
lp = -(nu+1)/2*log(1+(x-mu).^2/(s2*(nu-2))) + lZ;
dlp = (nu+1)*(mu-x)./(mu^2-2*mu*x+s2*(nu-2)+x.^2);