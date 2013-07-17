function K = covPeriodic(hyp, x, z, i)

% Stationary covariance function for a smooth periodic function, with period p:
%
% k(x,y) = sf2 * exp( -2*sin^2( pi*||x-y||/p )/ell^2 )
%
% where the hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)
%         log(sqrt(sf2)) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2011-01-05.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

n = size(x,1);
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