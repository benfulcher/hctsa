function K = covGaborard(hyp, x, z, i)

% Gabor covariance function with length scales ell and periods p. The covariance
% function is parameterized as:
%
% k(x,z) = h(x-z), h(t) = exp(-sum(t.^2./(2*ell.^2)))*cos(2*pi*sum(t./p)).
%
% The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D)
%         log(p_1)
%         log(p_2)
%          ..
%         log(p_D) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Note that covSM implements a weighted sum of Gabor covariance functions, but
% using an alternative (spectral) parameterization.
%
% Copyright (c) by Hannes Nickisch, 2014-09-26.
%
% See also COVFUNCTIONS.M, COVGABORISO, COVSM.M.

if nargin<2, K = '2*D'; return; end                        % report no of params
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);                                                % dimensionality
ell = exp(hyp(1:D));                                              % length scale
p = exp(hyp(D+1:2*D));                                                  % period

if dg                                              % compute squared distance d2
  d2 = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    d2 = sq_dist(diag(1./ell)*x');
  else                                                   % cross covariances Kxz
    d2 = sq_dist(diag(1./ell)*x',diag(1./ell)*z');
  end
end

dp = zeros(size(d2));                                % init sum(t)/p computation
if ~dg
  if xeqz                                                 % symmetric matrix Kxx
    for d=1:D, dp = dp + (x(:,d)*ones(1,size(x,1))-ones(n,1)*x(:,d)')/p(d); end
  else                                                   % cross covariances Kxz
    for d=1:D, dp = dp + (x(:,d)*ones(1,size(z,1))-ones(n,1)*z(:,d)')/p(d); end
  end
end

K = exp(-d2/2);
if nargin<4                                                         % covariance
  K = cos(2*pi*dp).*K;
else                                                               % derivatives
  d = mod(i-1,D)+1;
  if dg, K = K*0; end                                     % handle diagonal case
  if xeqz                                                 % symmetric matrix Kxx
    dd = x(:,d)*ones(1,size(x,1))-ones(n,1)*x(:,d)';
  else                                                   % cross covariances Kxz
    dd = x(:,d)*ones(1,size(z,1))-ones(n,1)*z(:,d)';
  end
  if i<=D                                                         % length scale
    K = (dd/ell(d)).^2.*cos(2*pi*dp).*K;
  elseif i<=2*D                                                         % period
    K = 2*pi*dd/p(d).*sin(2*pi*dp).*K;
  else
    error('Unknown hyperparameter')
  end
end