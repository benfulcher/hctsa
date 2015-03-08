function K = covPoly(d, hyp, x, z, i)

% Polynomial covariance function. The covariance function is parameterized as:
%
% k(x^p,x^q) = sf^2 * ( c + (x^p)'*(x^q) )^d 
%
% The hyperparameters are:
%
% hyp = [ log(c)
%         log(sf)  ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<3, K = '2'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

c = exp(hyp(1));                                          % inhomogeneous offset
sf2 = exp(2*hyp(2));                                           % signal variance
if d~=max(1,fix(d)), error('only nonzero integers allowed for d'), end  % degree

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

if nargin<5                                                        % covariances
  K = sf2*( c + K ).^d;
else                                                               % derivatives
  if i==1
    K = c*d*sf2*( c + K ).^(d-1);
  elseif i==2
    K = 2*sf2*( c + K ).^d;
  else
    error('Unknown hyperparameter')
  end
end