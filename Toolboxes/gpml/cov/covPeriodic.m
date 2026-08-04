function [K,dK] = covPeriodic(hyp, x, z)

% Stationary covariance function for a smooth periodic function, with period p 
% in 1d (see covPERiso and covPERard for multivariate data):
%
% k(x,z) = sf^2 * exp( -2*sin^2( pi*(x-z)/p )/ell^2 )
%
% where the hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-24.
%
% See also covFunctions.m.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D>1, error('Covariance is defined for 1d data only.'), end
ell = exp(hyp(1)); p = exp(hyp(2)); sf2 = exp(2*hyp(3));   % extract hyperparams

% precompute deviations and exploit symmetry of sin^2
if dg                                                               % vector txx
  T = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Txx
    T = pi/p*bsxfun(@plus,x,-x');
  else                                                   % cross covariances Txz
    T = pi/p*bsxfun(@plus,x,-z');
  end
end

S2 = (sin(T)/ell).^2;
K = sf2*exp( -2*S2 );                                              % covariances

if nargout>1
  dK = @(Q) dirder(Q,K,S2,T,ell,xeqz,x,z,p);      % directional hyper derivative
end

function [dhyp,dx] = dirder(Q,K,S2,T,ell,xeqz,x,z,p)
  Q = K.*Q; P = sin(2*T).*Q;
  dhyp = [4*(S2(:)'*Q(:)); 2/ell^2*(P(:)'*T(:)); 2*sum(Q(:))];
  if nargout > 1
    R = P./T; R(T==0) = 0;
    q2 = sum(R,2); q1 = sum(R,1)';
    if xeqz
      y = bsxfun(@times,q1+q2,x) - (R+R')*x;
    else
      Rz = R*z; y = bsxfun(@times,q2,x) - Rz;
    end
    dx = -2*pi^2/(ell*p)^2 * y;
  end