function [K,dK] = covOU(i, hyp, x, z)
% covOU i-times integrated Ornstein-Uhlenbeck process covariance function
%
%  [K,dK] = covOU(i, hyp, x, z)
%
% The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf)
%         log(sf0) ]
%
% For i=0, this considers the stochastic differential equation
%
%   ell * f'(x) + f(x) = N(0,sf^2), f(0) = N(f0,sf0^2), x>=0
%
% where 1/ell>0 is the decay rate and sf, sf0 are noise levels.
% N(m,v) is a Gaussian random variable with mean m and variance v.
%
% For i=0, this is the Ornstein-Uhlenbeck process covariance, for i=1, and
% assuming sf=sf0, this is the integrated Ornstein-Uhlenbeck process.
%
% The covariance function -- given that x,z>=0 -- is specified as:
% 
% i=0: k(x,z) = sf^2 * exp( -abs(x-z)/ell ) + (sf0^2-sf^2)*exp( -(x+z)/ell )
% i=1: k(x,z) = sf^2 * ell^2 * ( 2*min(x,z)/ell + r(x,z) ),
%              where r(x,z) = exp(-x/ell) + exp(-z/ell) - exp(-abs(x-z)/ell) - 1
%
% For more help on design of covariance functions, try "help covFunctions".
%
% References:
% See 10.2 of "Statistical Analysis of Stochastic Processes in Time",
% by J. K. Lindsey, CUP 2004.
% See "A Stochastic Model for Analysis of Longitudinal AIDS Data",
% by Taylor, Cumberland and Sy, JASA 1994.
%
% Copyright (c) by Juan Pablo Carbajal and Hannes Nickisch, 2017-01-19.
%
% See also cov/covW.m, covFunctions.m.

if nargin<3, K = '3'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D~=1, error('Covariance is defined for 1d data only.'), end
if any(x<0), error('Covariance is defined for nonnegative data only.'), end

ell = exp(hyp(1)); sf = exp(hyp(2)); sf0 = exp(hyp(3)); % obtain hyperparameters

ex = exp(-x/ell); ox = ones(size(x));                               % precompute
if dg                                                               % vector kxx
  M = x; G = 2-2*ex; dGa = -2*ex.*x/ell; E=[]; I=[]; ez=[]; oz=[]; 
  D = zeros(n,1); S = 2*x;
else
  if xeqz                                                 % symmetric matrix Kxx
    ez = ex; oz = ox; z = x;
  else                                                   % cross covariances Kxz
    ez = exp(-z/ell); oz = ones(size(z));
  end
  D = bsxfun(@minus,x,z'); S = bsxfun(@plus,x,z'); M = ox*z'; T = x*oz';
  I = bsxfun(@le,x,z'); M(I) = T(I); E = exp(-abs(D)/ell);
  G = 1+E-ex*oz'-ox*ez';
  dGa = (E.*abs(D)-(ex.*x)*oz'-ox*(ez.*z)')/ell;
end
if i==0
  eD = exp(-abs(D)/ell); eS = exp(-S/ell); K = sf^2*eD + (sf0^2-sf^2)*eS;
  dK = @(Q) dirder0(Q,D,S,eD,eS,ell,sf,sf0,dg,xeqz);
else
  K = 2*ell*sf^2*( M - ell/2*G );
  dK = @(Q) dirder1(Q,K,E,I,D,G,M,dGa,ex,ez,ox,oz,ell,sf,dg,xeqz);
end

function [dhyp,dx] = dirder0(Q,D,S,eD,eS,ell,sf,sf0,dg,xeqz)
  qes = Q(:)'*eS(:); qed = Q(:)'*eD(:);
  R = sf^2*eD.*abs(D) + (sf0^2-sf^2)*eS.*S;
  dhyp = [Q(:)'*R(:)/ell; 2*sf^2*(qed-qes); 2*sf0^2*qes];
  if nargout > 1
    if dg, dx = zeros(size(dpx));
    else   dx =   (sf^2-sf0^2)*sum(Q .*eS,2) - sf^2*sum(Q .*eD.*sign(D),2);
      if xeqz
        dx = dx + (sf^2-sf0^2)*sum(Q'.*eS,2) - sf^2*sum(Q'.*eD.*sign(D),2);
      end
    end
    dx = dx/ell;
  end

function [dhyp,dx] = dirder1(Q,K,E,I,D,G,M,dGa,ex,ez,ox,oz,a,sf,dg,xeqz)
  R = 2*a*M - a^2*(2*G+dGa);
  dhyp = [sf^2*(R(:)'*Q(:)); 2*(K(:)'*Q(:)); 0];
    if dg
      dx = zeros(size(ex));
    else
      dx = sum(Q.*I,2) + sum(Q.*(E.*sign(D)-ex*oz'),2)/2; % M+G
      if xeqz
        dx = dx + sum(Q.*(1-I),1)';                       % M
        dx = dx + sum(Q.*(E.*sign(-D)-ox*ez'),1)'/2;      % G
      end
      dx = 2*a*sf^2 * dx;
    end
