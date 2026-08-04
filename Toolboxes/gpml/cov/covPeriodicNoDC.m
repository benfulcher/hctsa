function [K,dK] = covPeriodicNoDC(hyp, x, z)

% Stationary covariance function for a smooth periodic function, with period p:
%
% k(x,z) = sf^2 * [k0(pi*(x-z)/p) - f(ell)] / [1 - f(ell)]
%        with k0(t) = exp( -2*sin^2(t)/ell^2 ) and f(ell) = \int 0..pi k0(t) dt.
%
% The constant (DC component) has been removed and marginal variance is sf^2.
% The hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)
%         log(sf) ]
%
% Note that covPeriodicNoDC converges to covCos as ell goes to infinity.
%
% Copyright (c) by James Robert Lloyd and Hannes Nickisch 2016-04-24.
%
% See also covFunctions.m, cov/covCos.m.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

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

K = covD(2*T,ell); K = sf2*K;                                       % covariance

if nargout>1
  ebi0 = embi0(1/ell^2); ebi1 = embi1(1/ell^2);
  dK = @(Q) dirder(Q,K,T,ell,p,sf2,xeqz,x,z,ebi0,ebi1);  % directional hyp deriv
end

function [dhyp,dx] = dirder(Q,K,T,ell,p,sf2,xeqz,x,z,ebi0,ebi1)
  S2 = (sin(T)/ell).^2;
  if ell>1e4                                              % limit for ell->infty
    Z = zeros(size(T));                    % no further progress in ell possible
  elseif 1/ell^2<3.75
    cK = cos(2*T); ecK = exp(cK/ell^2);
    b0 = besseli(0,1/ell^2);
    b1 = besseli(1,1/ell^2);
    Z =    2*(exp(1/ell^2)-ecK    )*b1 ...
         - 2*(exp(1/ell^2)-ecK.*cK)*b0 ...
         + 4*exp(2*(cos(T)/ell).^2).*sin(T).^2;
    Z = sf2/(ell*(exp(1/ell^2)-b0))^2 * Z;
  else
    cK = cos(2*T); ecK = exp((cK-1)/ell^2);
    b0 = ebi0; b1 = ebi1;
    Z =    2*(1-ecK)*b1 - 2*(1-ecK.*cK)*b0 ...
         + 4*exp(2*(cos(T).^2-1)/ell^2).*sin(T).^2;
    Z = sf2/(ell*(1-b0))^2 * Z;
  end
  if ell>1e4                                              % limit for ell->infty
    Y = 2*sf2*                         sin(2*T).*T;
    a = 1;
  elseif 1/ell^2<3.75
    c = 1/(exp(1/ell^2)-b0)/ell^2;
    Y = 2*c*sf2*exp( cos(2*T)   /ell^2).*sin(2*T).*T;
    a = 1/(1-b0*exp(-1/ell^2))/ell^2;
  else
    c = 1/(1-b0)/ell^2;
    Y = 2*c*sf2*exp((cos(2*T)-1)/ell^2).*sin(2*T).*T;
    a = c;
  end
  dhyp = [Z(:)'*Q(:); Y(:)'*Q(:); 2*(Q(:)'*K(:))];
  if nargout > 1
    Kdc = sf2*exp( -2*S2 ); % inkl. DC component
    R = Kdc.*sin(2*T).*Q./T; R(T==0) = 0;
    r2 = sum(R,2); r1 = sum(R,1)';
    if xeqz
      y = bsxfun(@times,r1+r2,x) - (R+R')*x;
    else
      Rz = R*z; y = bsxfun(@times,r2,x) - Rz;
    end
    dx = -2*a*pi^2/p^2 * y;
  end

function K = covD(D,ell)                   % evaluate covariances from distances
  if ell>1e4                                              % limit for ell->infty
    K = cos(D);
  elseif 1/ell^2<3.75
    K = exp(cos(D)/ell^2);
    b0 = besseli(0,1/ell^2);
    K = (K-b0)/(exp(1/ell^2)-b0);
  else
    K = exp((cos(D)-1)/ell^2);
    b0 = embi0(1/ell^2);
    K = (K-b0)/(1-b0);
  end

function y = embi0(x)      % = exp(-x)*besseli(0,x) => 9.8.2 Abramowitz & Stegun
  y = 3.75/x;
  y = 0.39894228     + 0.01328592*y   + 0.00225319*y^2 - 0.00157565*y^3 ...
    + 0.00916281*y^4 - 0.02057706*y^5 + 0.02635537*y^6 - 0.01647633*y^7 ...
    + 0.00392377*y^8;
  y = y/sqrt(x);

function y = embi1(x)      % = exp(-x)*besseli(1,x) => 9.8.4 Abramowitz & Stegun
  y = 3.75/x;
  y = 0.39894228     - 0.03988024*y   - 0.00362018*y^2 + 0.00163801*y^3 ...
    - 0.01031555*y^4 + 0.02282967*y^5 - 0.02895312*y^6 + 0.01787654*y^7 ...
    - 0.00420059*y^8;
  y = y/sqrt(x);