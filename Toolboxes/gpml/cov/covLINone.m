function K = covLINone(hyp, x, z, i)

% Linear covariance function with a single hyperparameter. The covariance
% function is parameterized as:
%
% k(x^p,x^q) = (x^p'*x^q + 1)/t^2;
%
% where the P matrix is t2 times the unit matrix. The second term plays the
% role of the bias. The hyperparameter is:
%
% hyp = [ log(t) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

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

if nargin<4                                                        % covariances
  K = it2*(1+K);
else                                                               % derivatives
  if i==1
    K = -2*it2*(1+K);
  else
    error('Unknown hyperparameter')
  end
end
