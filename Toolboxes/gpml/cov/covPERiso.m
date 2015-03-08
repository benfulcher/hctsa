function K = covPERiso(cov, hyp, x, z, i)

% Stationary periodic covariance function for an isotropic stationary covariance
% function k0 such as covMaterniso, covPPiso, covRQiso and covSEiso.
% Isotropic stationary means that the covariance function k0(x,z) depends on the
% data points x,z only through the squared distance
% dxz = (x-z)'*inv(P)*(x-z) where the P matrix is ell^2 times the unit matrix.
% The covariance function is parameterized as:
%
% k(x,z) = k0(u(x),u(z)), u(x) = [sin(pi*x/p); cos(pi*x/p)]
%
% where the period p belongs to covPERiso and hyp0 belong to k0:
%
% hyp = [ log(p)
%         hyp0 ]
%
% The first hyperparameter of k0 is the log lengthscale hyp0(1) = log(ell).
% Note that for k0 = covSEiso and D = 1, a faster alternative is covPeriodic.
%
% Copyright (c) by Hannes Nickisch, 2013-10-15.
%
% See also COVFUNCTIONS.M.

nocov = false;                   % default case when no cov argument is provided
if nargin==0, cov = {@covSEiso}; nocov = true; end                % default case
if isnumeric(cov)       % detect old version where the cov parameter was missing
  % i <- z, z <- x, x <- hyp, hyp <- cov
  if nargin>3, i = z; end
  if nargin>2, z = x; end
  if nargin>1, x = hyp; end
  hyp = cov; cov = {@covSEiso}; nocov = true;
end

if nocov && nargin<2 || ~nocov && nargin<3         % report number of parameters
  K = ['(1+',feval(cov{:}),')']; return
end
if nocov && nargin<3 || ~nocov && nargin<4, z = []; end    % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
p = exp(hyp(1));

if nocov && nargin<4 || ~nocov && nargin<5
  [x,z] = u(x,z,p,dg);                      % apply the embedding u:IR^D->IR^2*D
  K = feval(cov{:},hyp(2:end),x,z);
else
  if i==1
    if dg                                                   % compute distance d
      d = zeros([n,1,D]);
    else
      if xeqz                                             % symmetric matrix Kxx
        d = repmat(reshape(x,n,1,D),[1,n, 1])-repmat(reshape(x,1,n, D),[n,1,1]);
      else                                               % cross covariances Kxz
        nz = size(z,1);
        d = repmat(reshape(x,n,1,D),[1,nz,1])-repmat(reshape(z,1,nz,D),[n,1,1]);
      end
    end
    d = 2*pi*d/p; dD2_dlp = -2*sum(sin(d).*d,3);        % derivative dD2/dlog(p)
    [x,z] = u(x,z,p,dg);                    % apply the embedding u:IR^D->IR^2*D
    if dg                                            % compute squared distances
      D2 = zeros(n,1);
    else
      if xeqz, D2 = sq_dist(x'); else D2 = sq_dist(x',z'); end
    end
    % reconstruct derivative w.r.t. D2 from derivative w.r.t. log(ell)
    dK_dD2 = feval(cov{:},hyp(2:end),x,z,1)./(-2*D2); dK_dD2(D2<1e-12) = 0;
    K = dK_dD2.*dD2_dlp;                                      % apply chain rule
  else
    [x,z] = u(x,z,p,dg);                    % apply the embedding u:IR^D->IR^2*D
    K = feval(cov{:},hyp(2:end),x,z,i-1);
  end
end

function [x,z] = u(x,z,p,dg)                % apply the embedding u:IR^D->IR^2*D
  x = 2*pi*x/p; x = [sin(x), cos(x)];
  if numel(z)>0 && ~dg, z = 2*pi*z/p; z = [sin(z), cos(z)]; end