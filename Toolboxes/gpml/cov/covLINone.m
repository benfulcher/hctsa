function [K,dK] = covLINone(hyp, x, z)

% Linear covariance function with a single hyperparameter. The covariance
% function is parameterized as:
%
% k(x,z) = (x'*z + 1)/t^2;
%
% where the P matrix is t2 times the unit matrix. The second term plays the
% role of the bias. The hyperparameter is:
%
% hyp = [ log(t) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-23.
%
% See also covFunctions.m.

if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

it2 = exp(-2*hyp);                                                  % t2 inverse

% precompute inner products
if dg                                                               % vector kxx
  K = sum(x.*x,2);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = x*x';
  else                                                   % cross covariances Kxz
    K = x*z';
  end
end
K = it2*(K+1);                                                     % covariances
if nargout > 1, dK = @(Q) dirder(Q,K,x,z,it2,xeqz,dg); end      % dir derivative

function [dhyp,dx] = dirder(Q,K,x,z,it2,xeqz,dg)
  dhyp = -2*Q(:)'*K(:);
  if nargout>1
    if dg
      dx = zeros(size(x));
    else
      if xeqz
        dx = it2*(Q*x+Q'*x);
      else
        dx = it2*Q*z;
      end
    end
  end