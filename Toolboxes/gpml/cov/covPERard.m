function K = covPERard(cov, hyp, x, z, i)

% Stationary periodic covariance function for a stationary covariance function
% k0 such as covMaternard, covPPard, covRQard and covSEard.
% Stationary means that the covariance function k0(x,z) depends on the
% data points x,z only through the squared distance
% dxz = (x-z)'*inv(P)*(x-z) where the P matrix is diagonal with ARD parameters
% ell_1^2,...,ell_D^2, where D is the dimension of the input space.
% The covariance function is parameterized as:
%
% k(x,z) = k0(u(x),u(z)), u(x) = [sin(pi*x/p); cos(pi*x/p)]
%
% where the period p belongs to covPERiso and hyp0 belong to k0:
%
% hyp = [ log(p_1)
%         log(p_2)
%          .
%         log(p_D)
%         hyp0 ]
%
% The first D hyperparameters of k0 are the log lengthscales such that
% hyp0(i) = log(ell(i)) for i=1..D.
% Note that for k0 = covSEard and D = 1, a faster alternative is covPeriodic.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.
%
% See also COVFUNCTIONS.M.

nocov = false;                   % default case when no cov argument is provided
if nargin==0, cov = {@covSEard}; nocov = true; end                % default case
if isnumeric(cov)       % detect old version where the cov parameter was missing
  % i <- z, z <- x, x <- hyp, hyp <- cov
  if nargin>3, i = z; end
  if nargin>2, z = x; end
  if nargin>1, x = hyp; end
  hyp = cov; cov = {@covSEard}; nocov = true;
end

if nocov && nargin<2 || ~nocov && nargin<3         % report number of parameters
  K = ['(D+',feval(cov{:}),')']; return
end
if nocov && nargin<3 || ~nocov && nargin<4, z = []; end    % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
p = exp(hyp(1:D));
lell = hyp(D+(1:D));
hyp0 =  [lell; lell; hyp(2*D+1:end)];

if nocov && nargin<4 || ~nocov && nargin<5
  [x,z] = u(x,z,p,dg);                      % apply the embedding u:IR^D->IR^2*D
  K = feval(cov{:},hyp0,x,z);
else
  if i<=D
    if dg                                                   % compute distance d
      di = zeros([n,1]);
    else
      if xeqz                                             % symmetric matrix Kxx
        di = repmat(reshape(x(:,i),n, 1),[1, n])...
            -repmat(reshape(x(:,i),1, n),[n, 1]);
      else                                               % cross covariances Kxz
        nz = size(z,1);
        di = repmat(reshape(x(:,i),n, 1),[1,nz])...
            -repmat(reshape(z(:,i),1,nz),[n, 1]);
      end
    end
    di = 2*pi*di/p(i); dD2_dlpi = -2*sin(di).*di;   % derivative dD2i/dlog(p(i))
    [x,z] = u(x,z,p,dg);                    % apply the embedding u:IR^D->IR^2*D
    if dg                                            % compute squared distances
      D2 = zeros(n,1);
    else
      if xeqz
        D2 = sq_dist(x(:,[i,i+D])');
      else
        D2 = sq_dist(x(:,[i,i+D])',z(:,[i,i+D])');
      end
    end
    % reconstruct derivative w.r.t. D2i from derivative w.r.t. log(ell(i))
    dK_dD2 = feval(cov{:},hyp0,x,z,i) + feval(cov{:},hyp0,x,z,i+D);
    dK_dD2 = dK_dD2./(-2*D2); dK_dD2(D2<1e-12) = 0;
    K = dK_dD2.*dD2_dlpi;                                     % apply chain rule
  else
    [x,z] = u(x,z,p,dg);                    % apply the embedding u:IR^D->IR^2*D
    if i<=2*D
      K = feval(cov{:},hyp0,x,z,i-D) + feval(cov{:},hyp0,x,z,i);
    else
      K = feval(cov{:},hyp0,x,z,i);
    end
  end
end

function [x,z] = u(x,z,p,dg)                % apply the embedding u:IR^D->IR^2*D
  x = 2*pi*x*diag(1./p); x = [sin(x), cos(x)];
  if numel(z)>0 && ~dg, z = 2*pi*z*diag(1./p); z = [sin(z), cos(z)]; end