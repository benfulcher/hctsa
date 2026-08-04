function [K,dK] = covCos(hyp, x, z)

% Stationary covariance function for a sinusoid with period p in 1d:
%
% k(x,z) = sf^2*cos(2*pi*(x-z)/p)
%
% where the hyperparameters are:
%
% hyp = [ log(p)
%         log(sf) ]
%
% Note that covPeriodicNoDC converges to covCos as ell goes to infinity.
%
% Copyright (c) by James Robert Lloyd and Hannes Nickisch, 2016-11-05.
%
% See also covFunctions.m, cov/covPeriodicNoDC.m.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D~=1, error('Covariance is defined for 1d data only.'), end
p   = exp(hyp(1));
sf2 = exp(2*hyp(2));

% precompute deviations and exploit symmetry of cos
if dg                                                               % vector txx
  T = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Txx
    T = 2*pi/p*bsxfun(@plus,x,-x');
  else                                                   % cross covariances Txz
    T = 2*pi/p*bsxfun(@plus,x,-z');
  end
end

K = sf2*cos(T);                                                    % covariances
if nargout > 1
  dK = @(Q) dirder(Q,K,T,x,p,sf2,dg,xeqz);
end

function [dhyp,dx] = dirder(Q,K,T,x,p,sf2,dg,xeqz)
  dhyp = [sf2*(sin(T(:)).*T(:))'*Q(:); 2*Q(:)'*K(:)];
  if nargout > 1
    R = -sf2*pi/p * Q .* sin(T);
    if dg
      dx = zeros(size(x));
    else
      if xeqz
        dx = 2*(sum(R,2)-sum(R,1)');
      else
        dx = 2*sum(R,2);
      end
    end
  end