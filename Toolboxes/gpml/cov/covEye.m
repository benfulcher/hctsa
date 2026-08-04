function [K,dK] = covEye(hyp, x, z)

% Independent covariance function, i.e. "white noise", with unit variance.
% The covariance function is specified as:
%
% k(x^p,x^q) = \delta(p,q)
%
% \delta(p,q) is a Kronecker delta function which is 1 iff p=q and zero
% otherwise in mode 1).
% In cross covariance mode 2) two data points x_p and z_q are considered equal
% if their difference norm |x_p-z_q| is less than eps, the machine precision.
% The hyperparameters are:
%
% hyp = [ ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2016-04-18.
%
% See also covFunctions.m.

tol = eps;   % threshold on the norm when two vectors are considered to be equal
if nargin<2, K = '0'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
dg = strcmp(z,'diag');                                          % determine mode

n = size(x,1);

if dg                                                               % vector kxx
  K = ones(n,1);
else
  if isempty(z)                                           % symmetric matrix Kxx
    K = eye(n);
  else                                                   % cross covariances Kxz
    K = double(sq_dist(x',z')<tol*tol);
  end
end

if nargout > 1
  dK = @(Q) dirder(Q,x);                          % directional hyper derivative
end

function [dhyp,dx] = dirder(Q,x)
  dhyp = zeros(0,1); if nargout > 1, dx = zeros(size(x)); end