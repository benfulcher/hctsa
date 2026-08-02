function varargout = covMatern(mode, par, d, varargin)

% Matern covariance function with nu = d/2 and isotropic distance measure. For
% d=1 the function is also known as the exponential covariance function or the 
% Ornstein-Uhlenbeck covariance in 1d. The covariance function is:
%
%   k(x,z) = f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3, f(t)=1+t+t^2/3 for d=5 and
%      f(t)=1+t+2*t^2/5+t^3/15  for d=7.
%
% The covariance function can also be expressed for non-integer d.
%
%   k(x,z) = 2^(1-nu)/gamma(nu) * s^nu * besselk(nu,s), s = sqrt(2*nu)*r
%
% Note that for d->oo the covariance converges to the squared exponential.
%
% Here r is the Mahalanobis distance sqrt(maha(x,z)). The function takes a
% "mode" parameter, which specifies precisely the Mahalanobis distance used, see
% covMaha. The function returns either the number of hyperparameters (with less
% than 3 input arguments) or it returns a covariance matrix and (optionally) a
% derivative function.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-05-29.
%
% See also cov/covSE.m, cov/covMaha.m.

if nargin < 1, error('Mode cannot be empty.'); end            % no default value
if nargin < 2, par = []; end                                           % default
varargout = cell(max(1, nargout), 1);                  % allocate mem for output
if nargin < 5, varargout{1} = covMaha(mode,par); return, end

if any(d==[1,3,5,7])
  switch d                                              % df(t) = (f(t)-f'(t))/t
    case 1, f = @(t) 1;                    df = @(t) 1./t;
    case 3, f = @(t) 1+t;                  df = @(t) 1;
    case 5, f = @(t) 1+t.*(1+t/3);         df = @(t) (1+t)/3;
    case 7, f = @(t) 1+t.*(1+t.*(6+t)/15); df = @(t) (1+t+t.^2/3)/5; 
  end
            m = @(t,f) f(t).*exp(-t); dm = @(t,f) df(t).*exp(-t);
     k = @(d2)              m(sqrt(d*d2),f);
  if d==1
    dk = @(d2,k) set_zero( -dm(sqrt(  d2),f)/2, d2==0 );  % fix limit case d2->0
  else
    dk = @(d2,k)           -dm(sqrt(d*d2),f)*d/2;
  end
else                                                      % use non-integer mode
  k = @(d2) psi(d2,d/2); dk = @(d2,k) dpsi(d2,k,d/2);
end

[varargout{:}] = covMaha(mode, par, k, dk, varargin{:});

function A = set_zero(A,I), A(I) = 0;

function k = psi(d2,nu) % = 2^(1-nu)/gamma(nu) * r^nu*besselk(nu,r) book Eq 4.14
  c = (1-nu)*log(2)-gammaln(nu);     % evaluate 2^(1-nu)/gamma(nu) in log domain
  r = sqrt(2*nu*d2); k = exp(c+nu*log(r)).*besselk(nu,r);
  k(r<1e-7) = 1;                     % fix lim_r->0, see Abramowitz&Stegun 9.6.9
  i = isnan(k) | isinf(k);                             % detect strange behavior
  if any(i(:)), k(i) = exp(-d2(i)/2); end               % fix limit_nu->oo covSE

function dk = dpsi(d2,k,nu)
  r = sqrt(2*nu*d2);
  dk = -nu*0.5^nu/gamma(nu)*r.^(nu-2).*...
         (r.*besselk(nu-1,r)-2*nu*besselk(nu,r)+r.*besselk(nu+1,r));
  i = isnan(dk) | isinf(dk);                           % detect strange behavior
  if any(i(:)), dk(i) = (-1/2)*k(i); end                % fix limit_nu->oo covSE