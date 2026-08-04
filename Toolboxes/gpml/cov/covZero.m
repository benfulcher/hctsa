function [K,dK] = covZero(hyp, x, z)

% Constant (degenerate) covariance function, with zero variance.
% The covariance function is specified as:
%
% k(x,z) = 0
%
% hyp = [ ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2016-04-17.
%
% See also covFunctions.m.

if nargin<2, K = '0'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
dg = strcmp(z,'diag');                                          % determine mode

n = size(x,1);

if dg                                                               % vector kxx
  K = zeros(n,1);
else
  if isempty(z)                                           % symmetric matrix Kxx
    K = zeros(n,n);
  else                                                   % cross covariances Kxz
    K = zeros(n,size(z,1));
  end
end

if nargout > 1
  dK = @(Q) dirder(Q,K,x);                        % directional hyper derivative
end

function [dhyp,dx] = dirder(Q,K,x)
  dhyp = zeros(0,1); if nargout > 1, dx = zeros(size(x)); end
