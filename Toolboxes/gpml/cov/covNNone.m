function [K,dK] = covNNone(hyp, x, z)

% Neural network covariance function with a single parameter for the distance
% measure. The covariance function is parameterized as:
%
% k(x,z) = sf2 * asin(x'*P*z / sqrt[(1+x'*P*x)*(1+z'*P*z)])
%
% where the x and z vectors on the right hand side have an added extra bias
% entry with unit value. P is ell^-2 times the unit matrix and sf2 controls the
% signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sqrt(sf2) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-23.
%
% See also covFunctions.m.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

ell2 = exp(2*hyp(1));
sf2 = exp(2*hyp(2));

sx = 1 + sum(x.*x,2);
if dg                                                               % vector kxx
  A = sx./(sx+ell2); sz = sx;
else
  if xeqz                                                 % symmetric matrix Kxx
    S = 1 + x*x'; sz = sx;
    A = S./(sqrt(ell2+sx)*sqrt(ell2+sx)');
  else                                                   % cross covariances Kxz
    S = 1 + x*z'; sz = 1 + sum(z.*z,2);
    A = S./(sqrt(ell2+sx)*sqrt(ell2+sz)');
  end
end

K = sf2*asin(A);                                                   % covariances
if nargout > 1
  dK = @(Q) dirder(Q,K,A,sx,sz,ell2,sf2,x,z,dg,xeqz);     % dir hyper derivative
end

function [dhyp,dx] = dirder(Q,K,A,sx,sz,ell2,sf2,x,z,dg,xeqz)
  n = size(x,1);
  if dg
    V = A;
  else
    vx = sx./(ell2+sx);
    if xeqz
      V = repmat(vx/2,1,n) + repmat(vx'/2,n,1);
    else  
      vz = sz./(ell2+sz); nz = size(z,1);
      V = repmat(vx/2,1,nz) + repmat(vz'/2,n,1);
    end
  end
  P = Q./sqrt(1-A.*A);
  dhyp = [-2*sf2*sum(sum((A-A.*V).*P)); 2*Q(:)'*K(:)];
  if nargout > 1
    if dg
      dx = zeros(size(x));
    else
      W = P./(sqrt(ell2+sx)*sqrt(ell2+sz)'); ssx = sqrt(ell2+sx);
      if xeqz, W = W+W'; z = x; ssz = ssx; else ssz = sqrt(ell2+sz); end
      dx = sf2*(W*z - bsxfun(@times,x,((W.*A)*ssz)./ssx));
    end
  end