function [K,dK] = covFBM(hyp, x, z)
% covFBM Fractional Brownian motion covariance with Hurst index h from (0,1).
%
%  [K,dK] = covFBM(hyp, x, z)
%
% The hyperparameters are:
%
%   hyp = [ log(sf)
%          -log(1/h-1) ]
%
% For h=1/2, this is the Wiener covariance, for h>1/2, the increments are
% positively correlated and for h<1/2 the increments are negatively correlated.
%
% The covariance function -- given that x,z>=0 -- is specified as:
%
% k(x,z) = sf^2 / 2 * ( |x|^(2h) + |z|^(2h) - |x-z|^(2h) ), where
%
% sf^2 is the signal variance.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2017-01-27.
%
% See also cov/covW.m, covFunctions.m.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x); ox = ones(n,1);
if D~=1, error('Covariance is defined for 1d data only.'), end
if any(x<0) || (~xeqz && any(z<0))
  error('Covariance is defined for nonnegative data only.')
end
sf = exp(hyp(1)); h = 1/(1+exp(-hyp(2)));

if dg                                                               % vector kxx
  X = x; Z = x;
else
  if xeqz                                                 % symmetric matrix Kxx
    X = x*ox'; Z = X';
  else                                                   % cross covariances Kxz
    oz = ones(size(z,1),1);
    X = x*oz'; Z = ox*z';
  end
end

K = sf^2 * (abs(X).^(2*h) + abs(Z).^(2*h) - abs(X-Z).^(2*h))/2;     % signal var

if nargout > 1
  dK = @(Q) dirder(Q,K,X,Z,h,sf,dg,xeqz);                    % directional deriv
end

function [dhyp,dx] = dirder(Q,K,X,Z,h,sf,dg,xeqz)
  Ax = abs(X); Az = abs(Z); Ad = abs(X-Z);
  R = sf^2 * (  fixnan(log(Ax).*Ax.^(2*h)) + fixnan(log(Az).*Az.^(2*h)) ...
              - fixnan(log(Ad).*Ad.^(2*h)) );
  dhyp = [2*(Q(:)'*K(:)); (Q(:)'*R(:))*(h-h^2)];               % hyperparameters
  if nargout > 1
    if dg
      dx = zeros(size(x));
    else
      if xeqz
        R = sign(X).*Ax.^(2*h-1) - sign(X-Z).*Ad.^(2*h-1);
        dx = 2*h*sum(R.*(Q+Q'),2);
      else
        R = sign(X).*Ax.^(2*h-1) - sign(X-Z).*Ad.^(2*h-1);
        dx = 2*h*sum(R.*Q,2);
      end
    end
    dx = sf^2*dx/2;                                            % signal variance
  end

function B = fixnan(A), B = A; B(isnan(B)) = 0;               % replace NaN by 0
