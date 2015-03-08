function K = covPeriodic(hyp, x, z, i)

% Stationary covariance function for a smooth periodic function, with period p 
% in 1d (see covPERiso and covPERard for multivariate data):
%
% k(x,z) = sf^2 * exp( -2*sin^2( pi*||x-z||/p )/ell^2 )
%
% where the hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2011-01-05.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D>1, error('Covariance is defined for 1d data only.'), end
ell = exp(hyp(1));
p   = exp(hyp(2));
sf2 = exp(2*hyp(3));

% precompute distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sqrt(sq_dist(x'));
  else                                                   % cross covariances Kxz
    K = sqrt(sq_dist(x',z'));
  end
end

K = pi*K/p;
if nargin<4                                                        % covariances
    K = sin(K)/ell; K = K.*K; K =   sf2*exp(-2*K);
else                                                               % derivatives
  if i==1
    K = sin(K)/ell; K = K.*K; K = 4*sf2*exp(-2*K).*K;
  elseif i==2
    R = sin(K)/ell; K = 4*sf2/ell*exp(-2*R.*R).*R.*cos(K).*K;
  elseif i==3
    K = sin(K)/ell; K = K.*K; K = 2*sf2*exp(-2*K);
  else
    error('Unknown hyperparameter')
  end
end