function [K,dK,D2] = covMaha(mode, par, k, dk, hyp, x, z)
% COVMAHA Mahalanobis distance-based covariance function.
%
% Report number of hyperparameters
%  s = COVMAHA (mode)
%  s = COVMAHA (mode, par)
%  s = COVMAHA (mode, par, k)
%  s = COVMAHA (mode, par, k, dk)
%  s = COVMAHA (mode, par, k, dk, hyp)
%
% Evaluation of k(x,x)
%  k = COVMAHA (mode, par, k, dk, hyp, x)
%  k = COVMAHA (..., x, [])
%
% Evaluation of diag(k(x,x))
%  k = COVMAHA (mode, par, k, dk, hyp, x, 'diag')
%
% Evaluation of k(x,z)
%  k = COVMAHA (mode, par, k, dk, hyp, x, z)
%
% Evaluation of function derivatives dk (w.r.t. hyp and x)
%  [k, dk] = COVMAHA (mode, par, k, dk, hyp, x)
%  [k, dk] = COVMAHA (..., x, [])
%  [k, dk] = COVMAHA (..., x, 'diag')
%  [k, dk] = COVMAHA (..., x, z)
%
% Call covFunctions to get an explanation of outputs in each mode.
%
% A third output argument in the evaluation mode returns the squared distance
%  [k, dk, r^2] = COVMAHA (..., x, z)
%
% The covariance function is parameterized as:
%
% k(x,z) = k(r^2)
%
% with
%
% r^2 = maha(x,P,z) = (x-z)' * inv(P) * (x-z),
%
% where the matrix P is the metric.
%
% Parameters:
% 1) mode, par:
% We offer different modes (mode) with their respective parameters (par):
% mode =   par =   inv(P) =         hyp =
%   'eye'    []      eye(D)           []
%   'iso'    []      ell^2*eye(D)     [log(ell)]
%   'ard'    []      diag(ell.^2)     [log(ell_1); ..; log(ell_D)]
%   'proj'   d       L'*L             [L_11; L_21; ..; L_dD]
%   'fact'   d       L'*L + diag(f)   [L_11; L_21; ..; L_dD; log(f_1); ..; log(f_D)]
%   'vlen'   llen    l(x,z)^2*eye(D)  [hyp_llen]
%
% In the last mode, the covariance function is turned into a nonstationary
% covariance by a variable lengthscale
%
% l(x,z) = sqrt(( len(x)^2 + len(z)^2 )/2),
%
% where len(x) = exp(llen(x)) and llen is provided as additional parameter
% in the form of a GPML mean function cell array. The final expression for the
% 'vlen' covariance is:
%
%   k(x,z) = ( len(x)*len(z)/l(x,z)^2 )^(D/2) * k( (x-z)'*(x-z)/l(x,z)^2 )
%
% 2) k, dk:
% The functional form of the covariance is governed by two functions:
%
% k:  r^2         --> k(x,z), r^2 = maha(x,P,z) = (x-z)' * inv(P) * (x-z)
% dk: r^2, k(x,z) --> d k(x,z) / d r2
%
% For example, the squared exponential covariance uses
%
%   k = @(r2) exp(-r2/2); dk = @(r2,k) (-1/2)*k;
%
% Note that not all functions k, dk give rise to a valid i.e. positive
% semidefinite covariance function k(x,z).
%
% 3) hyp, x, z:
% These input parameters follow the usual covariance function interface. For the
% composition of hyp, see 1).
%
% 4) K, dK:
% See the usual covariance function interface.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% See also COVFUNCTIONS

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2017-01-11.

if nargin<1, mode = 'eye'; end, if nargin <2, par = []; end     % default values
mode_list = '''eye'', ''iso'', ''ard'', ''proj'', ''fact'', or ''vlen''';
switch mode
  case 'eye',  ne = '0';
  case 'iso',  ne = '1';
  case 'ard',  ne = 'D';
  case 'proj', ne = [num2str(par),'*D'];
  case 'fact', ne = [num2str(par),'*D+D'];
  case 'vlen', ne = feval(par{:});
  otherwise,   error('Parameter mode is either %s.',mode_list)
end

if nargin<6, K = ne; return; end                   % report number of parameters
if nargin<7, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');             % sort out different modes
[n,D] = size(x); ne = eval(ne);                                     % dimensions
hyp = hyp(:);                                 % make sure hyp is a column vector
if numel(hyp)~=ne, error('Wrong number of hyperparameters'), end

switch mode                                           % mvm with metric A=inv(P)
  case {'eye','vlen'}
               A = @(x) x;
               dAdhyp = @(dAdiag,dAmvm) zeros(0,1);
  case 'iso',  A = @(x) x*exp(-2*hyp);
               dAdhyp = @(dAdiag,dAmvm) -2*sum(A(dAdiag'));
  case 'ard',  A = @(x) bsxfun(@times,x,exp(-2*hyp'));
               dAdhyp = @(dAdiag,dAmvm) -2*A(dAdiag')';
  case 'proj', d = par; L = reshape(hyp,d,D); A = @(x) (x*L')*L;
               dAdhyp = @(dAdiag,dAmvm) 2*reshape(dAmvm(L')',d*D,1);
  case 'fact', d = par; L = reshape(hyp(1:d*D),d,D); f = exp(hyp(d*D+1:end));
               A = @(x) (x*L')*L + bsxfun(@times,x,f');
               dAdhyp = @(dAdiag,dAmvm)[2*reshape(dAmvm(L')',d*D,1); f.*dAdiag];
end

[D2,dmaha] = maha(x,A,z); T = 1; L2 = 1;          % compute Mahalanobis distance
if isequal(mode,'vlen')                         % evaluate variable lengthscales
  lx = exp(feval(par{:},hyp,x));                 % L2 = (lx^2+lz^2)/2, P = lx*lz
  if dg
    L2 = lx.*lx; P = L2;
  else
    if xeqz, lz = lx; else lz = exp(feval(par{:},hyp,z)); end
    L2 = bsxfun(@plus,(lx.*lx)/2,(lz.*lz)'/2); P = lx*lz';
  end
  D2 = D2./L2; T = (P./L2).^(D/2);                   % non-stationary covariance
end
K = k(D2);                                                 % evaluate covariance
if nargout > 1                                          % directional derivative
  dK = @(Q) dirder(Q,K,dk,T,D2,L2,dmaha,dAdhyp,mode,par,hyp,x,z,dg,xeqz);
end
if isequal(mode,'vlen'), K = K.*T; end

function [dhyp,dx] = dirder(Q,K,dk,T,D2,L2,dmaha,dAdhyp, mode,par,hyp,x,z,dg,xeqz)
  R = T.*dk(D2,K).*Q;
  switch mode
    case 'eye',  dx = dmaha(R); dhyp = zeros(0,1);               % fast shortcut
    case 'iso',  dx = dmaha(R); dhyp = -2*R(:)'*D2(:);           % fast shortcut
    case 'vlen', dx = dmaha(R); lx2 = exp(-2*feval(par{:},hyp,x));
      dx = bsxfun(@times,dx,lx2); D = size(x,2);  % only correct for lx = const.
      if dg
        dhyp = zeros(size(hyp));
      else
        [llx,dllx] = feval(par{:},hyp,x); lx = exp(llx);
        if xeqz, lz = lx; dllz = dllx;
        else
          [llz,dllz] = feval(par{:},hyp,z); lz = exp(llz);
        end
        A=(D/2)*Q.*K.*((lx*lz')./L2).^(D/2-1)./L2; B=(D2.*R+A.*(lx*lz'))./L2;
        dhyp = dllx(lx.*(A*lz-sum(B,2).*lx)) + dllz(lz.*(A'*lx-sum(B,1)'.*lz));
      end
    otherwise, [dx,dAdiag,dAmvm] = dmaha(R); dhyp = dAdhyp(dAdiag,dAmvm);
  end

% Mahalanobis squared distance function for A spd
% D2 = maha(x,A,z) = (x-z)'*A*(x-z)
% dx(Q) = d tr(Q'*D2) / d x
% dA(Q) = d tr(Q'*D2) / d A
function [D2,dmaha] = maha(x,A,z)
  if nargin<2, A = @(x) x; end                     % make sure the metric exists
  if nargin<3, z = []; end                                  % make sure z exists
  xeqz = isempty(z); dg = strcmp(z,'diag');           % sort out different modes
  n = size(x,1); m = size(z,1);                                     % dimensions
  if dg                                                            % vector d2xx
    D2 = zeros(n,1);
  else
    % Computation of a^2 - 2*a*b + b^2 is less stable than (a-b)^2 because
    % numerical precision can be lost when both a and b have very large absolute
    % value and the same sign. For that reason, we subtract the mean from the
    % data beforehand to stabilise the computations. This is OK because the
    % squared error is independent of the mean.
    if xeqz
      mu = mean(x,1);
    else
      mu = (m/(n+m))*mean(z,1) + (n/(n+m))*mean(x,1); z = bsxfun(@minus,z,mu);
    end
    x = bsxfun(@minus,x,mu);
    Ax = A(x); sax = sum(x.*Ax,2);
    if xeqz                                              % symmetric matrix D2xx
      Az = Ax; saz = sax;
    else                                                      % cross terms D2xz
      Az = A(z); saz = sum(z.*Az,2);
    end                % remove numerical noise at the end and ensure that D2>=0
    D2 = max(bsxfun(@plus,sax,bsxfun(@minus,saz',2*x*Az')),0);     % computation
  end
  if nargout>1, dmaha = @(Q) maha_dirder(Q,x,A,z); end

function [dx,dAdiag,dAmvm] = maha_dirder(Q,x,A,z)       % directional derivative
  if nargin<3, z = []; end
  xeqz = isempty(z); dg = strcmp(z,'diag');           % sort out different modes
  q2 = sum(Q,2); q1 = sum(Q,1)'; sym = @(X) (X+X')/2;
  dAdense = size(x,2)<5;    % estimated break-even between O(D*D*n) and O(4*D*n)
  if dg
    dx = zeros(size(x));
    if nargout > 1
      dAdiag = zeros(size(x,2),1);
      dAmvm = @(r) zeros(size(r));
    end
  else
    if xeqz
      y = bsxfun(@times,q1+q2,x) - (Q+Q')*x;
      if nargout > 1
        if dAdense         % construct a dense matrix dA of size DxD in O(D*D*n)
          dA = sym(x'*y); dAdiag = diag(dA); dAmvm = @(r) dA*r;
        else          % just perform an MVM avoiding a DxD matrix dA in O(4*D*n)
          dAdiag = sum(x.*y,1)'; dAmvm = @(r) (x'*(y*r) + y'*(x*r))/2;
        end
      end
    else
      Qz = Q*z; y = bsxfun(@times,q2,x) - Qz;
      if nargout > 1
        yz = y-Qz; qz = bsxfun(@times,q1,z);
        if dAdense         % construct a dense matrix dA of size DxD in O(D*D*n)
          dA = sym(x'*yz + z'*qz);
          dAdiag = diag(dA); dAmvm = @(r) dA*r;
        else             % just perform an MVM avoiding a DxD matrix in O(8*D*n)
          dAdiag = sum(x.*yz,1)'+sum(z.*qz,1)';
          dAmvm = @(r) (x'*(yz*r)+yz'*(x*r) + z'*(qz*r)+qz'*(z*r))/2;
        end
      end
    end
    dx = 2*A(y);
  end
