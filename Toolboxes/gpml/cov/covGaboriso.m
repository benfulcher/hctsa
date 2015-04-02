function K = covGaboriso(hyp, x, z, i)

% Gabor covariance function with length scale ell and period p. The 
% covariance function is parameterized as:
%
% k(x,z) = h(x-z) with h(t) = exp(-t'*t/(2*ell^2))*cos(2*pi*sum(t)/p).
%
% The hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)   ]
%
% Note that covSM implements a weighted sum of Gabor covariance functions, but
% using an alternative (spectral) parameterization.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2014-09-26.
%
% See also COVFUNCTIONS.M, COVGABORARD.M, COVSM.M.

if nargin<2, K = '2'; return; end                          % report no of params
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);                                                % dimensionality
ell = exp(hyp(1));                                                % length scale
p = exp(hyp(2));                                                        % period

if dg                                              % compute squared distance d2
  d2 = zeros(n,1);
else
  if xeqz                                                 % symmetric matrix Kxx
    d2 = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    d2 = sq_dist(x'/ell,z'/ell);
  end
end

dp = zeros(size(d2));                                % init sum(t)/p computation
if ~dg
  if xeqz                                                 % symmetric matrix Kxx
    for d=1:D, dp = dp + (x(:,d)*ones(1,size(x,1))-ones(n,1)*x(:,d)')/p; end
  else                                                   % cross covariances Kxz
    for d=1:D, dp = dp + (x(:,d)*ones(1,size(z,1))-ones(n,1)*z(:,d)')/p; end
  end
end

K = exp(-d2/2);
if nargin<4                                                        % covariances
  K = cos(2*pi*dp).*K;
else                                                               % derivatives
  if i==1                                                         % length scale
    K = d2.*cos(2*pi*dp).*K;
  elseif i==2                                                           % period
    K = 2*pi*dp.*sin(2*pi*dp).*K;
  else
    error('Unknown hyperparameter')
  end
end