function [lp,dlp] = priorInvGauss(mu,lam,x)

% Univariate Inverse Gaussian hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(-lam*(x-mu)^2/(2*mu^2*x)) / sqrt(2*pi*x^3/lam)
%
% where mu(1x1) is the mean parameter, lam(1x1) is the scale parameter
% and x(1xN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also PRIORDISTRIBUTIONS.M.

if nargin<2, error('mu and lam parameters need to be provided'), end
if ~(isscalar(mu)&&isscalar(lam))
  error('mu and lam parameters need to be scalars'),end
if nargin<3                                                    % return a sample
  n = randn; y = n*n;
  r = mu + mu*(mu*y-sqrt(4*mu*lam*y+mu^2*y^2))/(2*lam);
  z = rand;
  if z<=mu/(mu+r)
    lp = r;
  else
    lp = mu^2/r;
  end
  return
end

lp  = -lam*(x-mu).^2./(2*mu^2*x) - log(2*pi*x.^3/lam)/2;
q = (x-mu)./x;
dlp = -lam*q.*(2-q)/(2*mu^2) - 3./(2*x);
lp(x<0) = -inf; dlp(x<0) = 0;