function [lp,dlp] = priorSmoothBox1(a,b,eta,x)

% Univariate smoothed box prior distribution with linear decay in the log domain
% and infinite support over the whole real axis.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%  p(x) = sigmoid(eta*(x-a))*(1-sigmoid(eta*(x-b))),
%         where sigmoid(z) = 1/(1+exp(-z))
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

if nargin<4 % inverse sampling
  u = exp((b-a)*eta*rand());
  lp = log((u-1)/(exp(-eta*a)-u*exp(-eta*b)))/eta;
  return
end

[lpa,dlpa] = logr(eta*(x-a)); [lpb,dlpb] = logr(-eta*(x-b));
 lp = lpa + lpb - log(b-a) + log(1-exp((a-b)*eta));
dlp = eta*(dlpa - dlpb);

% r(z) = 1/(1+exp(-z)), log(r(z)) = -log(1+exp(-z))
function [lr,dlr] = logr(z)
  lr = z; ok = -35<z;  lr(ok) = -log(1+exp(-z(ok)));
 dlr = 1./(1+exp(z));