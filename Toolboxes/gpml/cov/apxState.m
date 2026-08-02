function varargout = apxState(cov,varargin)

% apxSparse - Covariance function for state space representations.
%
% This is a special covariance function which handles some non-standard
% requirements of sparse covariances
% 
% [1] Solin, Sarkka, Explicit link between periodic covariance functions and
%     state space models, AISTATS, 2014.
% [2] Sarkka, Solin, Hartikainen, Spatiotemporal learning via infinite-dim
%     Bayesian filtering and smoothing. IEEE Signal Processing Magazine, 2013.
% [3] Hartikainen, Sarkka, Kalman filtering and smoothing solutions to temporal
%     GP regression models, MLSP, 2010.
%
% Adapted from code by Arno Solin in the GPstuff toolbox
%                             http://research.cs.aalto.fi/pml/software/gpstuff/.
%
% Copyright (c) by Hannes Nickisch, 2018-06-01.
%
% See also covFunctions.m, cov/apx.m, inf/infLaplace.m, inf/infGaussLik.m.

derivtest = false;

if nargin<1, error('Not enough parameters provided.'), end
varargout = cell(nargout, 1);        % allocate right number of return arguments

% mode  1) perform Kalman filtering
if     strcmp(cov,'kalman')
  [varargout{:}] = kalman(varargin{:}); return
% mode  2) perform Kalman smoothing
elseif strcmp(cov,'rts')
  [varargout{:}] = rts(varargin{:}); return
% mode  3) compute transition matrices
elseif strcmp(cov,'trans')
  [varargout{:}] = trans(varargin{:});  return
end

if nargin>2, [varargout{:}] = feval(cov{:},varargin{:}); return, end % pure eval
cov = f2s(cov);   % function handles -> strings+'_' to avoid extern dependencies
if ischar(cov) || isa(cov,'function_handle'), cov  = {cov};  end     % make cell
if nargin<2, varargout{1} = feval(cov{:}); return, end   % number of hyperparams 
[varargout{:}] = feval(cov{:},varargin{:});           % set up state space model

if derivtest                              % switch to perform derivative testing
  hyp = varargin{1};
  dev = @(x,z) norm(x(:)-z(:))/max([1,norm(x(:)),norm(z(:))]);
  nh = numel(hyp); h = 1e-5;
  [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = feval(cov{:},hyp);
  dFn = zeros(size(dF)); dQcn = zeros(size(dQc)); dPinfn = zeros(size(dPinf));
  for i=1:nh
    hyph = hyp; hyph(i) = hyph(i)+h;
    [Fi,L,Qci,H,Pinfi] = feval(cov{:},hyph);
    dFn(:,:,i) = (Fi-F)/h; dQcn(:,:,i) = (Qci-Qc)/h;
    dPinfn(:,:,i) = (Pinfi-Pinf)/h;
  end
  err = [dev(dF,dFn), dev(dQc,dQcn), dev(dPinf,dPinfn)];
  assert(err(1)<6*h)
  assert(err(2)<7*h)
  assert(err(3)<6*h)
end

% ldB2 = sum(log(z))/2
% dz = d z / d theta, where dAs = d As / d theta etc.
% t(f) = exp(tnu*f-ttau/2*f^2) ttau = 1/s2, tnu = mu/s2 with natural params
%
% allow for assumed density filtering (ADF) using moment matching
% in that case tnu,ttau are updated
function [z,ga, tnu,ttau, ms,Ps,dz] = kalman(As,Qs, H,Pinf, tnu,ttau, t,tid, ...
                                                             mom, dAs,dQs,dPinf)
  d = numel(H); m = zeros(d,1); P = Pinf; n = size(tnu,1); nt = numel(t); % init
  rm = nargout>4; if rm, ms = zeros(d,nt);   end
  rp = nargout>5; if rp, Ps = zeros(d,d,nt); end
  deriv = nargout>6;                                 % Are derivatives required?
  if deriv, nh = size(dPinf,3); dz = zeros(n,nh); dP = dPinf; end % allocate mem
  deriv = deriv && nargin>11;                        % Are derivatives possible?  
  adf = nargin>8 && numel(mom)>0; % Are we performing assumed density filtering?
  z = zeros(n,1); sym = @(X) X+X'; ga = zeros(n,1);
  for j=1:nt
    l = tid(j);
    A = As(:,:,j); Q = Qs(:,:,j);              % use the method by Davison & Man
    if j>1 && deriv
      for i=1:nh
        dP(:,:,i) = sym(dAs(:,:,i,j)*P*A') + A*dP(:,:,i)*A' + dQs(:,:,i,j);
      end
    end              % shift from f_t-1 towards f_t using f_t ~ N(A_t*f_t-1,Q_t)
    if j>1, m = A*m; P = A*P*A'+Q; end % prediction = propagate mean m and cov P
    if l<=n                   % Kalman filter update for training instances only
      fmu = H*m; w = P*H'; fs2 = H*w;     % latent marginal, cavity distribution
      if adf                                                       % perform ADF
        [lZ,dlZ,d2lZ] = mom(fmu,fs2,l);        % prop moments through likelihood
        ttau(l) = -d2lZ/(1+d2lZ*fs2);                  % perform moment matching
        tnu(l)  = (dlZ-fmu*d2lZ)/(1+d2lZ*fs2);
        ttau(l) = max(ttau(l),0); % enforce positivity->lower bound ttau by zero
      end
      z(l) = ttau(l)*fs2+1; % Gauss: H*m~N(r,W^-1) ttau = W, tnu = W*r, r=y-m(x)
      k = ttau(l)*w/z(l);                        P = P-k*w';  % variance-related
      c = ttau(l)*fmu - tnu(l); ga(l) = -c/z(l); m = m+ga(l)*w;   % mean-related
%       % beta(l) = -c/sqrt(ttau(l)*z(l)), where beta = chol(L)'\r      
      if deriv
        for i=1:nh
          dw = dP(:,:,i)*H'; dz(l,i) = ttau(l)*H*dw; 
          dk = dw/z(l) - w*dz(l,i)/z(l)^2;
          dP(:,:,i) = dP(:,:,i) - ttau(l)*dk*w' - k*dw';
        end
      end
    end
    if rm, ms(:,j) = m; end, if rp, Ps(:,:,j) = P; end
  end

function [ms,Ps,rho] = rts(As,Qs,ms,Ps,H) % Rauch-Tung-Striebel Kalman smoothing
  [d,nt] = size(ms); m = ms(:,nt); P = Ps(:,:,nt); rho = zeros(nt,1);     % init
  for i=nt:-1:2                                        % run backward iterations
    A = As(:,:,i); Q = Qs(:,:,i); mi = ms(:,i-1); Pi = Ps(:,:,i-1);
    Pj = A*Pi*A'+Q; [L,npd] = chol(Pj,'lower');  % Cholesky otherwise add jitter
    if npd>0
      e = eig(Pj); s = 1e-5*max(e); 
      C = (Pi*A')/(Pj+s*diag(1+rand(d,1)));
    else
      C = Pi*A'/L'/L;
    end
    dmi = C*(m-A*mi); m = mi+dmi; P = Pi+C*(P-Pj)*C';                   % update
    ms(:,i-1) = m; Ps(:,:,i-1) = P;                             % store smoothed
    if nargin>4 && nargout>2, rho(i-1) = H*dmi; end             % only if needed
  end

function [As,Qs, t,tid, dAs,dQs] = trans(F,L,Qc,H,Pinf,stable,x,ng,dF,dQc,dPinf)
  deriv = nargin>8 && nargout>4;       % decide whether derivatives are required
  x0 = min(x); if stable, x0 = x0-0.1; else x0 = 0; end  % compute virtual point
  d = numel(H); [t,tid] = sort([x;x0]); nt = numel(t);   % include virtual point
  if nargin<8 || numel(ng)==0, ng = nt; end                  % set default value
  ng = min(ng,nt);      % do not interpolate if exact computation costs the same
  As = zeros(d,d,ng); Qs = zeros(d,d,ng); dt = diff(t(:));     % allocate memory
  if deriv, nh = size(dF,3); dAs = zeros(d,d,nh,ng); dQs = zeros(d,d,nh,ng); end
  if ng==nt                                                        % set up grid
    dg = [0;dt];                                                    % exact case
  else
    dg = linspace(min(dt),max(dt),ng-4); step = dg(2)-dg(1);   % equispaced case
    dg = [dg(1)-[2,1]*step, dg, dg(end)+[1,2]*step]';
  end
  dg_new = inf; A = zeros(d); Q = zeros(d);                        % init memory
  for i=1:ng
    dg_old = dg_new; dg_new = dg(i);               % g is sorted and hence dg>=0
    if abs(dg_new-dg_old)>1e-9               % only recompute if spacing changes
      if deriv
        [A,Q, dA,dQ] = integ(F,L,Qc,dg(i),Pinf,stable, dF,dQc,dPinf);
      else
        [A,Q] = integ(F,L,Qc,dg(i),Pinf,stable);
      end
    end
    if deriv, dAs(:,:,:,i) = dA; dQs(:,:,:,i) = dQ; end
    As(:,:,i) = A; Qs(:,:,i) = Q;
  end
  if ng~=nt
    W = apxGrid('interp',{dg},[0;dt],3); W(1,:) = 0;  % cubic conv interpolation
    if deriv
      dAs = reshape(reshape(dAs,d^2*nh,ng)*W',d,d,nh,nt);
      dQs = reshape(reshape(dQs,d^2*nh,ng)*W',d,d,nh,nt);
    end
    As = reshape(reshape(As,d^2,ng)*W',d,d,nt);
    Qs = reshape(reshape(Qs,d^2,ng)*W',d,d,nt);
  end
  if deriv, dQs(:,:,:,1) = dPinf; end                       % set initial values
  As(:,:,1) = eye(d); Qs(:,:,1) = Pinf;

function [A,Q, dA,dQ] = integ(F,L,Qc,dt,Pinf,stat, dF,dQc,dPinf)
  deriv = nargin>8 && nargout>2;       % decide whether derivatives are required
  d = size(F,1);
  if deriv
    Z = zeros(d); nh = size(dF,3); dA = zeros(d,d,nh); dQ = zeros(d,d,nh);
    for i=1:nh
      % see eq 10+11 in Najfeld and Havel, Derivatives of the matrix exponential
      % and their computation, AAM 16, p 321-375, 1995.
      if stat          % chain rule through matrix exponentiation of double size
        E = expm(dt*[F,Z; dF(:,:,i),F]);
        A = E(1:d,1:d); dAi = E(d+1:2*d,1:d); dA(:,:,i) = dAi;
        Q = Pinf - A*Pinf*A'; S = A*Pinf*dAi'; S = S+S';
        dQ(:,:,i) = dPinf(:,:,i) - A*dPinf(:,:,i)*A' - S;
      else            % closed-form integration by matrix fraction decomposition
        W = L*dQc(:,:,i)*L'; U = L*Qc*L'; G = dF(:,:,i);
        E = expm(dt*[F U Z Z; Z -F' Z Z; G W F U; Z -G' Z -F']);
        A = E(1:d,1:d); dA(:,:,i) = E(2*d+1:3*d,1:d);
        E1 = E(    1:  d,d+1:2*d); E2 = E(  d+1:2*d,d+1:2*d);
        E3 = E(2*d+1:3*d,d+1:2*d); E4 = E(3*d+1:4*d,d+1:2*d);
        Q = E1/E2; dQ(:,:,i) = (E3 - Q*E4)/E2;
      end
    end
  else
    % solve Q = int_[t=0..dt] P(dt-t)*L*Qc*L'*P(dt-t)', where P(t) = expm(t*F)
    if stat                                          % see Davison and Man, 1968
      A = expm(dt*F); Q = Pinf - A*Pinf*A';   % discrete-time state trans matrix
    else              % closed-form integration by matrix fraction decomposition
      E = expm(dt*[F,L*Qc*L'; zeros(d),-F']);                     % for LTI SDEs
      A = E(1:d,1:d); Q = E(1:d,d+1:2*d)/E(d+1:2*d,d+1:2*d);
    end
  end

function covs = f2s(covf)   % replace function handles by strings in nested cell
  if isa(covf,'function_handle')
    covs = func2str(covf);
    if strncmp('cov',covs,3), covs = [covs,'_']; end
  elseif ischar(covf) && strncmp('cov',covf,3)
    covs = [covf,'_'];
  elseif iscell(covf)
    covs = cell(size(covf));
    for i=1:numel(covf), covs{i} = f2s(covf{i}); end
  else
    covs = covf;
  end

function varargout = covSum_(varargin)              % wrapper, allocate and call
  varargout = cell(nargout,1); [varargout{:}] = covSumProd(false,varargin{:});

function varargout = covProd_(varargin)             % wrapper, allocate and call
  varargout = cell(nargout,1); [varargout{:}] = covSumProd(true,varargin{:});

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covSumProd(isprod,cov,hyp)
  nc = numel(cov);                      % number of terms in covariance function
  for ii = 1:nc                              % iterate over covariance functions
    f = cov(ii); if iscell(f{:}), f = f{:}; end % expand cell array if necessary
    j(ii) = cellstr(feval(f{:}));                        % collect number hypers
  end
  if nargin<3                                      % report number of parameters
    F = char(j(1)); for ii=2:nc, F = [F, '+', char(j(ii))]; end, return
  end
  v = [];             % v vector indicates to which covariance parameters belong
  for ii = 1:nc, v = [v repmat(ii, 1, eval(char(j(ii))))]; end, stat = true;
  if isprod                                % different init for product/sum mode
    F = 0;  dF = []; L = 1;  Qc = 1;  dQc = 1;  H = 1;  Pinf = 1;  dPinf = 1;
  else
    F = []; dF = []; L = []; Qc = []; dQc = []; H = []; Pinf = []; dPinf = [];
  end
  for ii = 1:nc                               % iteration over summand functions
    f = cov(ii); if iscell(f{:}), f = f{:}; end % expand cell array if necessary
    [Fi,Li,Qci,Hi,Pinfi,stati, dFi,dQci,dPinfi] = ...
                                   feval(f{:}, hyp(v==ii)); stat = stat & stati;
    if isprod                                                     % product mode
      [F1,dF1] = dkron(F,dF, eye(size(Fi)),zeros(size(dFi)));
      [F2,dF2] = dkron(eye(size(F)),zeros(size(dF)), Fi,dFi);
      F = F1+F2; dF = dF1+dF2;
      L = kron(L,Li);
      [Qc,dQc] = dkron(Qc,dQc,Qci,dQci);
      H = kron(H,Hi);
      [Pinf,dPinf] = dkron(Pinf,dPinf,Pinfi,dPinfi);
    else                                                              % sum mode
      F = blkdiag(F,Fi);
      L = blkdiag(L,Li); Qc   = blkdiag(Qc,Qci);
      H = [H,Hi];        Pinf = blkdiag(Pinf,Pinfi);
      dF = blkdiag3(dF,dFi); dQc = blkdiag3(dQc,dQci);
      dPinf = blkdiag3(dPinf,dPinfi);
    end
  end

function [C,dC] = dkron(A,dA,B,dB)
  if numel(A)==1 && A==1, C = B; dC = dB; return, end             % trivial case
  C = kron(A,B);
  nha = size(dA,3)*(numel(dA)>0); nhb = size(dB,3)*(numel(dB)>0);
  dC = zeros(size(C,1),size(C,2),nha+nhb);
  for i=1:nha, dC(:,:,i) = kron(dA(:,:,i),B); end
  for i=1:nhb, dC(:,:,nha+i) = kron(A,dB(:,:,i)); end

function C = blkdiag3(A,B)
  sa = [size(A,1),size(A,2),size(A,3)]; if numel(A)==0, sa(3) = 0; end
  sb = [size(B,1),size(B,2),size(B,3)]; if numel(B)==0, sb(3) = 0; end
  sc = sa+sb;
  C = zeros(sc); if ~isempty(A), C(1:sa(1),1:sa(2),1:sa(3))=A; end  % dimensions
  C(sa(1)+1:end,sa(2)+1:end,sa(3)+1:end) = B;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covLINiso_(hyp)
  if nargin<1, F = '1'; return, end           % return number of hyperparameters
  Qc = 0; sf = exp(-hyp(1));
  F = [0 1;0 0]; dF = zeros(2); L = [0;1]; dQc = 0; H = [1,0];
  global t0, if numel(t0)==0, t0 = 0; end % ugly but the only way to get this in
  t0 = 0;
  Pinf = sf^2*[t0^2,t0; t0,1]; dPinf = -2*Pinf; stat = false;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covNoise_(hyp)
  if nargin<1, F = '1'; return, end           % return number of hyperparameters
  ell = 1e-9; sf = exp(hyp(1)); Qc = 2; a = 1;    % lim_[ell->0] of covMaterniso
  [F,dF,L,Qc,dQc,H, Pinf,dPinf] = covAtomic(Qc,a,ell,sf); stat = true;
  dF = dF(:,:,2); dQc = dQc(:,:,2); dPinf = dPinf(:,:,2);

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covConst_(hyp)
  if nargin<1, F = '1'; return, end           % return number of hyperparameters
  sf = exp(hyp(1)); stat = true;
  F = 0; dF = 0; L = sf; Qc = 0; dQc = 0; H = 1;
  Pinf = sf^2; dPinf = 2*sf^2;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covOU_(o,hyp)
  if nargin<2, F = '3'; return, end           % return number of hyperparameters
  a = 1; Qc = 2;
  ell = exp(hyp(1)); sf = exp(hyp(2));
  [F,dF,L,Qc,dQc,H] = covAtomic(Qc,a,ell,sf);
  Pinf = exp(2*hyp(3)); stat = false;
  dPinf = zeros(1,1,3); dPinf(3) = 2*Pinf; dF(:,:,3) = 0; dQc(:,:,3) = 0;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covW_(o,hyp)
  if nargin<2, F = '1'; return, end           % return number of hyperparameters
  if o<0, error('Order needs to be a non-negative integer.'), end
  a = zeros(1,o+1);
  Qc = 1; ell = 1; sf = exp(hyp(1));
  [F,dF,L,Qc,dQc,H] = covAtomic(Qc,a,ell,sf); dF = dF(:,:,2); dQc = dQc(:,:,2);
  Pinf = zeros(o+1); dPinf = zeros(o+1); stat = false;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covMaterniso_(o,hyp)
  if nargin<2, F = '2'; return, end           % return number of hyperparameters
  d = (o+1)/2;
  a = [1,cumprod((d-(1:d-1)+1)./(1:d-1))];    % binomial coefficients of (x+1)^d
  a = a .* (sqrt(o).^(d:-1:1));
  Qc = 2*sqrt(pi)*gamma((1+o)/2)/gamma(o/2)*sqrt(o)^o;   % spectral density at 0
  ell = exp(hyp(1)); sf = exp(hyp(2));
  [F,dF,L,Qc,dQc,H, Pinf,dPinf] = covAtomic(Qc,a,ell,sf); stat = true;

% squared exponentional scale mixture approximation of degree M
function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covRQiso_(hyp)
  M = 7;                                                  % approximation degree
  if nargin<1, F = '3'; return, end           % return number of hyperparameters
  ell = exp(hyp(1)); sf = exp(hyp(2)); alpha = exp(hyp(3));% get hyperparameters
  [x,w] = gaulag(M,alpha-1);                               % Gauss-Laguerre rule
  a = sqrt(w/gamma(alpha)); b = sqrt(alpha./x);           % scale mixture params
  cov = cell(M,1); for i=1:M, cov{i} = 'covSEiso_'; end
  hyp = [log(b*ell)'; log(a*sf)'];
  [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covSum_(cov,hyp);
  [da2,db] = dal(M,alpha);               % derivatives of a^2 and b w.r.t. alpha
  dla = alpha*da2./(2*a.^2); dlb = alpha*(db./b);
  dF = csm(dF,dla,dlb); dQc = csm(dQc,dla,dlb); dPinf = csm(dPinf,dla,dlb);

function Y = csm(X,dla,dlb)             % alpha chain rule through scale mixture
  sx = size(X,1); nx = size(X,3);
  Y = zeros(sx,sx,3);
  Y(:,:,1) = sum(X(:,:,1:2:nx),3);
  Y(:,:,2) = sum(X(:,:,2:2:nx),3);
  Y(:,:,3) = reshape(reshape(X(:,:,1:2:nx),[],nx/2)*dlb,sx,sx) + ...
             reshape(reshape(X(:,:,2:2:nx),[],nx/2)*dla,sx,sx);

function [x,w] = gaulag(M,beta)  % points/weights for generalized Gauss-Laguerre
  p = lag(M,beta); x = roots(p); p = lag(M+1,beta);
  w = gamma(M + beta + 1)*x ./ (factorial(M)*(M+1)^2*polyval(p,x).^2);

function c = lag(M,beta)          % generalized Laguerre polynomial coefficients
  i = 1:M; a = (2*i-1) + beta; b = sqrt(i(1:M-1) .* ((1:M-1) + beta));
  c = (-1)^M/factorial(M) * poly(diag(a)+diag(b,1)+diag(b,-1));

function [da2,db] = dal(M,alpha)      % derivatives of a^2 and b w.r.t. to alpha
  bin = @(n,k) gamma(n+1)./gamma(k+1)./gamma(n-k+1);      % binomial coefficient
  Ln = @(n,alpha) (-1).^(n:-1:0).* ...                    % generalised Laguerre
                          bin(n+alpha,0:n)./factorial(n:-1:0);
  dLn = @(n,alpha) (-1).^(n:-1:0).* ...                    % derivative of above
                          bin(n+alpha,0:n)./factorial(n:-1:0).* ...
                          (digamma(alpha+n+1) - digamma(alpha+1+(n:-1:0)));
  c  = Ln(M, alpha-1); x = roots(c);
  dx = polyval(dLn(M, alpha-1),x)./ polyval(Ln(M-1, alpha),x);
  c2 = Ln(M+1, alpha-1); dc2 = dLn(M+1, alpha-1);
  da2 = (gamma(M+alpha)*(polyval(c2,x).*dx ...
       + x.*(polyval(c2,x).*(digamma(alpha+M) - digamma(alpha)) ...
       + 2*polyval(Ln(M,alpha),x).*dx - 2*polyval(dc2,x))))./ ...
       ((1 + M)^2*factorial(M)*gamma(alpha).*polyval(c2,x).^3);
  db = (alpha./x).^(3/2)/(2*alpha^2).*(x-alpha*dx);

% Taylor approximation of degree N
function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covSEiso_(hyp)
  N = 7;                                                  % approximation degree
  if nargin<1, F = '2'; return, end           % return number of hyperparameters
  fn = factorial(N); Qc = sqrt(2*pi)*fn*2^N;
  p = zeros(1,2*N+1); % polynomial coeffs
  for n=0:N, p(end-2*n) = fn*2^(N-n)/factorial(n)/(-1)^(n); end
  r = roots(p); a = poly(r(real(r)<0)); a = a(end:-1:2);
  ell = exp(hyp(1)); sf = exp(hyp(2));                 % extract hyperparameters
  [F,dF,L,Qc,dQc,H, Pinf,dPinf] = covAtomic(Qc,a,ell,sf); stat = true;

function [F,L,Qc,H,Pinf,stat, dF,dQc,dPinf] = covPeriodic_(hyp)
  if nargin<1, F = '3'; return, end           % return number of hyperparameters
  N = 6;                                                  % approximation degree
  ell = exp(hyp(1)); p = exp(hyp(2)); sf2 = exp(2*hyp(3)); % extract hyperparams
  b = besseli(0:N,1/ell^2); a = 2*sf2*exp(-ell^-2)*b; a(1) = .5*a(1);   % coeffs
  da = 2*sf2/ell^2/exp(ell^-2)*[b(1)-b(2),2*((1+(1:N)*ell^2).*b(2:N+1)-b(1:N))];
  F = kron(diag(0:N),[0 -2*pi/p; 2*pi/p 0]); L = eye(2*(N+1));
  Qc = zeros(2*(N+1)); H = kron(ones(1,N+1),[1 0]);
  Pinf = kron(diag(a),eye(2)); stat = true;
  dF = zeros(2*(N+1),2*(N+1),3);
  dF(:,:,2) = -F;
  dQc = zeros(2*(N+1),2*(N+1),3);
  dPinf = zeros(2*(N+1),2*(N+1),3);
  dPinf(:,:,1) = kron(diag(da),eye(2));
  dPinf(:,:,3) = 2*Pinf;

function [F,dF,L,Qc,dQc,H, Pinf,dPinf] = covAtomic(Qc,a,ell,sf)
  % note: the inputs Qc and a are supposed to not depend on ell, sf
  if nargin<4, sf = 1; end, if nargin<3, ell = 1; end             % default case
  d = numel(a); a = a./ell.^(d:-1:1);                       % standard embedding
  F = diag(ones(d-1,1),1); F(d,:) = -a;                        % feedback matrix
  dF = zeros(d,d,2); dF(d,:,1) = a.*(d:-1:1);
  Qc = Qc/ell^(2*(d-1)+1);             % scale input  i.e. spectral density at 0
  Qc = sf^2*Qc;                        % scale output i.e. spectral density at 0
  dQc = cat(3,-(2*(d-1)+1)*Qc,2*Qc);
  H = [1,zeros(1,d-1)]; L = [zeros(d-1,1);1];  % observation, noise effect model
  if nargout>6  % Pinf=stationary covariance matrix solving the Riccati equation
    vec = @(X) X(:);    % solve Lyapunov equation: F*Pinf + Pinf*F' + L*Q*L' = 0
    M = kron(F,eye(d))+kron(eye(d),F); r = vec(L*Qc*L'); % solution costs O(d^6)
    if norm(r)<eps, Pinf = zeros(d,d); else Pinf = -reshape(M\r,d,d); end
    dPinf = zeros(d,d,2);
    for i=1:2                    % derivative via the implicit function therorem
      dM = kron(dF(:,:,i),eye(d))+kron(eye(d),dF(:,:,i));
      dr = vec(L*dQc(:,:,i)*L') + dM*vec(Pinf);
      if norm(dr)<eps, dPinf(:,:,i) = zeros(d,d); else
        dPinf(:,:,i) = -reshape(M\dr,d,d);
      end
    end  
  end