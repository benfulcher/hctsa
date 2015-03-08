function K = covLINiso(hyp, x, z, i)

% Linear covariance function with Automatic Relevance Determination (ARD). The
% covariance function is parameterized as:
%
% k(x^p,x^q) = x^p'*inv(P)*x^q
%
% where the P matrix is ell^2 times the unit matrix. The hyperparameters are:
%
% hyp = [ log(ell) ]
%
% Note that there is no bias term; use covConst to add a bias.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-10-13.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

ell = exp(hyp(1));
x = x/ell;

% precompute inner products
if dg                                                               % vector kxx
  K = sum(x.*x,2);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = x*x';
  else                                                   % cross covariances Kxz
    z = z/ell;
    K = x*z';
  end
end

if nargin>3                                                        % derivatives
  if i==1
    if dg
      K = -2*sum(x.*x,2);
    else
      if xeqz
        K = -2*x*x';
      else
        K = -2*x*z';
      end
    end
  else
    error('Unknown hyperparameter')
  end
end