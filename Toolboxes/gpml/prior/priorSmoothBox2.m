function [lp,dlp] = priorSmoothBox2(a,b,eta,x)

% Univariate smoothed box prior distribution with quadratic decay in the log
% domain and infinite support over the whole real axis.
% The plateau pdf is built by cutting a Gaussian into two parts and
% inserting a uniform distribution from a to b.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%                            / N(x|a,s^2),  x<=a,
%  p(x) =  1/(w*(1/eta+1)) * | 1         ,  a<x<b
%                            \ N(x|b,s^2),  b<=x,
%    where s = w/(eta*sqrt(2*pi)), w = abs(b-a)
%
% a(1x1) is the lower bound parameter, b(1x1) is the upper bound parameter,
% eta(1x1)>0 is the slope parameter and  x(1xN) contains query hyperparameters
% for prior evaluation. Larger values of eta make the distribution more
% box-like.
%
%          /------------\
%         /              \
% -------- |            | --------> x
%          a            b
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Jose Vallet and Hannes Nickisch, 2014-09-08.
%
% See also PRIORDISTRIBUTIONS.M.

if nargin<3, error('a, b and eta parameters need to be provided'), end
if b<=a, error('b must be greater than a.'), end
if ~(isscalar(a)&&isscalar(b)&&isscalar(eta))
  error('a, b and eta parameters need to be scalar values')
end

w = abs(b-a); sab = w/(eta*sqrt(2*pi));           % box width and boundary slope

if nargin<4                                                    % return a sample
  if rand<eta/(eta+1)
    lp = a+w*rand;
  else
    g = sab*randn; if g<0, lp = g+a; else lp = g+b; end
  end
  return
end

lp = zeros(size(x)) - log(w*(1/eta+1));
dlp = zeros(size(x));
i = x<a;
lp(i) = lp(i) - (x(i)-a).^2/(2*sab^2);
dlp(i) = (a-x(i))/sab^2;
i = x>b;
lp(i) = lp(i) - (x(i)-b).^2/(2*sab^2);
dlp(i) = (b-x(i))/sab^2;