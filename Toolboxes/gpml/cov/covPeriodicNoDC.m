function K = covPeriodicNoDC(hyp, x, z, i)

% Stationary covariance function for a smooth periodic function, with period p:
%
% k(x,x') = sf^2 * [k0(pi*(x-x')/p) - f(ell)] / [1 - f(ell)]
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
% Copyright (c) by James Robert Lloyd and Hannes Nickisch 2013-10-21.
%
% See also COVFUNCTIONS.M, COVCOS.M.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

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

K = 2*pi*K/p;
if nargin<4                                                        % covariances
  K = sf2*covD(K,ell);
else
  if i==1
    if ell>1e4                                            % limit for ell->infty
      K = zeros(size(K));                  % no further progress in ell possible
    elseif 1/ell^2<3.75
      cK = cos(K); ecK = exp(cK/ell^2);
      b0 = besseli(0,1/ell^2);
      b1 = besseli(1,1/ell^2);
      K =    2*(exp(1/ell^2)-ecK    )*b1 ...
           - 2*(exp(1/ell^2)-ecK.*cK)*b0 ...
           + 4*exp(2*(cos(K/2)/ell).^2).*sin(K/2).^2;
      K = sf2/(ell*(exp(1/ell^2)-b0))^2 * K;
    else
      cK = cos(K); ecK = exp((cK-1)/ell^2);
      b0 = embi0(1/ell^2);
      b1 = embi1(1/ell^2);
      K =    2*(1-ecK)*b1 - 2*(1-ecK.*cK)*b0 ...
           + 4*exp(2*(cos(K/2).^2-1)/ell^2).*sin(K/2).^2;
      K = sf2/(ell*(1-b0))^2 * K;
    end
  elseif i==2
    if ell>1e4                                            % limit for ell->infty
      K = sf2*sin(K).*K;
    elseif 1/ell^2<3.75
      K = exp(cos(K)/ell^2).*sin(K).*K;
      K = sf2/ell^2/(exp(1/ell^2)-besseli(0,1/ell^2))*K;
    else
      K = exp((cos(K)-1)/ell^2).*sin(K).*K;
      K = sf2/ell^2/(1-embi0(1/ell^2))*K;
    end
  elseif i==3
      K = 2*sf2*covD(K,ell);
  else
    error('Unknown hyperparameter')
  end
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