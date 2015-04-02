function K = covLINard(hyp, x, z, i)

% Linear covariance function with Automatic Relevance Determination (ARD). The
% covariance function is parameterized as:
%
% k(x^p,x^q) = x^p'*inv(P)*x^q
%
% where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
% D is the dimension of the input space. The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D) ]
%
% Note that there is no bias term; use covConst to add a bias.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = 'D'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

ell = exp(hyp);
[n,D] = size(x);
x = x*diag(1./ell);

% precompute inner products
if dg                                                               % vector kxx
  K = sum(x.*x,2);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = x*x';
  else                                                   % cross covariances Kxz
    z = z*diag(1./ell);
    K = x*z';
  end
end

if nargin>3                                                        % derivatives
  if i<=D
    if dg
      K = -2*x(:,i).*x(:,i);
    else
      if xeqz
        K = -2*x(:,i)*x(:,i)';
      else
        K = -2*x(:,i)*z(:,i)';
      end
    end
  else
    error('Unknown hyperparameter')
  end
end