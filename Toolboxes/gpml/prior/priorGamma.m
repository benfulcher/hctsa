function [lp,dlp] = priorGamma(k,t,x)

% Univariate Gamma hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(x/t)/gamma(k)*x^(k-1)/t^k
%
% where k(1x1) is the shape parameter, t(1x1) is the scale parameter
% and x(1xN) contains query hyperparameters for prior evaluation.
%
% Sampling is done using the algorithm from p. 53 of the paper Generating gamma
% variates by a modified rejection technique by J.H. Ahrens and U. Dieter, 1982,
% ACM, 25, 47–54.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also priorDistributions.m.

if nargin<2, error('k and t parameters need to be provided'), end
if ~(isscalar(k)&&isscalar(t))
  error('k and t parameters need to be scalars'), end
if nargin<3
  m = 1;
  d = k-floor(k); % fractional part
  v0 = exp(1)/(exp(1)+d);
  while true
    v = rand(3,1);
    if v(1)<=v0
      r = v(2)^(1/d); s = v(3)*r^(d-1);
    else
      r = 1-log(v(2)); s = v(3)*exp(-r);
    end
    if s<=r^(d-1)*exp(-r), break, end
    m = m+1;
  end
  u = rand(floor(k),1);
  lp = t*(r-sum(log(u)));
  return
end

lx = log(x);
lp  = -gammaln(k) - k*log(t) + (k-1)*lx - x/t;
dlp = (k-1)./x - 1/t;
lp(x<0) = -inf; dlp(x<0) = 0;