function K = covPref(cov, hyp, x, z, varargin)

% covPref - covariance function for preference learning. The covariance
% function corresponds to a prior on f(x1) - f(x2).
%
% k(x,z) = k_0(x1,z1) + k_0(x2,z2) - k_0(x1,z2) - k_0(x2,z1).
%
% The hyperparameters are:
%
% hyp = [ hyp_k0 ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% See Collaborative Gaussian Processes for Preference Learning, NIPS 2014.
%
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2014-11-01.
%
% See also COVFUNCTIONS.M.

if nargin<3, K = strrep(feval(cov{:}),'D','D/2'); return; end     % no of params

x1 = x(:,1:end/2); x2 = x(:,1+end/2:end);
if nargin<4 || isempty(z), z = x; end
if strcmp(z,'diag')
  n = size(x,1); diag12 = zeros(n,1);
  for i=1:n, diag12(i) = feval(cov{:},hyp,x1(i,:),x2(i,:),varargin{:}); end
  K =   feval(cov{:},hyp,x1,z,varargin{:}) ...
      + feval(cov{:},hyp,x2,z,varargin{:}) - 2*diag12;
else
  z1 = z(:,1:end/2); z2 = z(:,1+end/2:end);
  K =   feval(cov{:},hyp,x1,z1,varargin{:}) ...
      + feval(cov{:},hyp,x2,z2,varargin{:}) ...
      - feval(cov{:},hyp,x1,z2,varargin{:}) ...
      - feval(cov{:},hyp,x2,z1,varargin{:});
end
