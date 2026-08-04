function [K,dK] = covULL(hyp, x, z)
% covULL Stationary covariance function for underdamped linear Langevin process
%
%  [K,dK] = covULL(hyp, x, z)
%
% where the hyperparameters are:
%
%  hyp = [ log(mu)
%          log(omega)
%          log(a) ]
%
% The covariance is obtained by filtering white noise through an underdamped 
% 2nd order system 
%
% m * f''(x) + c * f'(x) + k * f(x) = N(0,sf^2).
%
% k(t) = a^2*exp(-mu*t) * ( sin(omega*t)/omega + cos(omega*t)/mu ),
%
% where t = abs(x-z), and
%
% mu = c/(2*m), omega = sqrt(k/m-mu^2), a = sf/(2*sqrt(m*k))
%
% Reference:
% See https://en.wikipedia.org/wiki/Harmonic_oscillator
%
% Copyright (c) by Robert MacKay, 2016-11-15.
%
% See also covFunctions.m.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D~=1, error('Covariance is defined for 1d data only.'), end

mu   = exp(hyp(1));
omega = exp(hyp(2));
a2 = exp(2*hyp(3));

% precompute deviations 
if dg                                                               % vector txx
  T = zeros(n,1);
else
  if xeqz                                                 % symmetric matrix Txx
    T = bsxfun(@minus,x,x');
  else                                                   % cross covariances Txz
    T = bsxfun(@minus,x,z');
  end
end

K = a2*exp(-mu*abs(T)) .* ( sin(omega*abs(T))/omega + cos(omega*T)/mu );   % cov
if nargout > 1
  dK = @(Q) dirder(Q,K,T,mu,omega,a2,dg,xeqz);
end

function [dhyp,dx] = dirder(Q,K,T,mu,omega,a2,dg,xeqz)
  A = mu*abs(T).*K + a2/mu*exp(-mu*abs(T)) .* cos(omega*T);
  B = abs(T).*cos(omega*T) - sin(omega*abs(T))/omega - omega*T.*sin(omega*T)/mu;
  B = a2*exp(-mu*abs(T)) .* B;
  dhyp = [-A(:)'*Q(:); B(:)'*Q(:); 2*(K(:)'*Q(:))];
  if nargout > 1
    R = -a2*(mu/omega+omega/mu)*exp(-mu*abs(T)) .* sin(omega*T) .* Q;
    if dg
      dx = zeros(size(Q,1),1);
    else
      if xeqz
        dx = sum(R,2)-sum(R,1)';
      else
        dx = sum(R,2);
      end
    end
  end
