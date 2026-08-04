function K = apx(hyp,cov,x,opt)

% (Approximate) linear algebra operations involving the covariance matrix K.
%
% A) Exact covariance computations.
%    There are no parameters in this mode.
%    Depending on the sign of W, we switch between
%    - the symmetric Cholesky mode [1], where B = I + sqrt(W)*K*sqrt(W), and
%    - the asymmetric LU mode [2],      where B = I + K*W.
%    Note that (1) is faster but it is only applicable if all W are positive.
%
% B) Sparse covariance approximations aka FITC [4], VFE [5] and SPEP [3].
%    We have a parameter opt.s from [0,1], the power in sparse power EP [3] 
%    interpolating between the Fully Independent Training Conditionals (FITC)
%    approach [4] and a Variational Free Energy (VFE) approximation [5].
%    In particular:
%  opt.s, default is 1.0 for FITC, opt.s = 0.0 corresponds to VFE.
%    Please see cov/apxSparse.m for details.
%
% C) Grid-based covariance approximations aka KISS-GP [6].
%    Please see cov/apxGrid.m for further details and more parameters.
%  opt.cg_tol,   default is 1e-6      as in Matlab's pcg function
%  opt.cg_maxit, default is min(n,20) as in Matlab's pcg function
%    The conjugate gradient-based linear system solver has two adjustable
%    parameters, the relative residual threshold for convergence opt.cg_tol and
%    the maximum number of MVMs opt.cg_maxit until the process stops.
%  opt.deg,      default is 3         degree of interpolation polynomial
%    For equispaced axes, opt.deg sets the degree of the interpolation
%    polynomial along each of the p axes. Here 0 means nearest neighbor,
%    1 means linear interpolation, and 3 uses a cubic.
%    For non-equispaced axes, only linear interpolation with inverse distance
%    weighting is offered and opt.deg is ignored.
%
%  opt.ldB2_method string indicating the possible modes (details below)
%      - 'scale' => (i)   scaled eigenvalue approach followed by Fiedler bound
%      - 'cheby' => (ii)  stochastic estimation using Chebyshev polynomials
%      - 'lancz' => (iii) stochastic estimation using Lanczos iterations
%  (i)   'scale'   the default method simply uses the eigenvalues of the
%         covariance matrix of the complete grid and rescales the values by the
%         ratio of number grid points N and number data points n; in case of
%         non-Gaussian inference i.e. W not being isotropic, we apply and
%         additional bounding step (Fiedler)
%  (ii)  'cheby'  employs Monte-Carlo trace estimation aka the Hutchinson method
%         and Chebyshev polynomials to approximate the term log(det(B))/2
%         stochastically, see [7]. The following parameters configure different
%         aspects of the estimator and are only valid if opt.ldB2_method='cheby'
%    opt.ldB2_hutch,        default is 10, number of samples for the trace estim
%    opt.ldB2_cheby_degree, default is 15, degree of Chebyshev approx polynomial
%    opt.ldB2_maxit,        default is 50, max # of MVMs to estimate max(eig(B))
%    opt.ldB2_seed,         default is [], random seed for the stoch trace estim
%  (iii) 'lancz' employs Monte-Carlo trace estimation aka the Hutchinson method
%         and Lanczos iterations with full Gram-Schmidt to approximate the
%         term log(det(B))/2 stochastically.
%        'lancz-arpack' is the same as above only that the ARPACK reverse
%         communication interface is used and partial orthogonalisation is used.
%         The following parameters configure different aspects of the estimator
%         and are only valid if opt.ldB2_method = 'lancz*'
%    opt.ldB2_hutch,        default is 10, number of samples for the trace estim
%    opt.ldB2_maxit,        default is 50, max # of MVMs per Lanczos run
%    opt.ldB2_seed,         default is [], random seed for the stoch trace estim
%
%  opt.stat = true returns a little bit of output summarising the exploited
%    structure of the covariance of the grid.
%    The log-determinant approximation employs Fiedler's 1971 inequality and a
%    rescaled version of the eigenspectrum of the covariance evaluated on the
%    complete grid.
%
% D) State space covariance approximations [8,9]
%    opt.deriv      flag allows to switch off derivative computations
%    opt.ngrid      specifies the number of grid points for matrix interpolation
%    opt.balance    flag to improve state space model numerical stability
%    opt.qridge     ridge added to Q matrices to improve condition, default 1e-7
%    opt.rridge     ridge added to R matrix to improve condition,   default 0
%    opt.alpha_alg  alorithm for solving linear systems,        default 'spingp'
%                   'spingp' : SpInGP algorithm using sparse matrix algebra
%                   'kalman' : Kalman filtering/smoothing
%
% The call K = apx(hyp,cov,x,opt) yields a structure K with a variety of
% fields.
% 1) Matrix-vector multiplication with covariance matrix
%    K.mvm(x) = K*x
% 2) Projection and its transpose (unity except for mode B) Sparse approx.)
%    post.alpha = K.P(solveKiW(f-m))
%    post.L = L = @(r) -K.P(solveKiW(K.Pt(r)))
% 3) Linear algebra functions depending on W
%    [ldB2,solveKiW,dW,dldB2,L,triB] = K.fun(W)
%   a) Log-determinant (approximation), called f in the sequel
%      ldB2 = log(det(B))/2
%   b) Solution of linear systems
%      solveKiW(r) = (K+inv(W)) \ r
%   c) Log-determinant (approximation) derivative w.r.t. W
%      dW = d f / d W, where f = ldB2(W), exact value dW = diag(inv(B)*K)/2
%   d) Log-determinant (approximation) derivative w.r.t. hyperparameters
%      dhyp = dldB2(alpha,a,b)
%      Q = d f / d K, exact value would be Q = inv(K+inv(W))
%      R = alpha*alpha' + 2*a*b'
%      Here dhyp(i) = tr( (Q-R)'*dKi )/2, where dKi = d K / d hyp(i).
%   e) Matrix (operator) to compute the predictive variance
%      L = -K.P(solveKiW(K.Pt(r))) either as a dense matrix or function.
%      See gp.m for details on post.L.
%   f) triB = trace(inv(B))
%
% [1] Seeger, GPs for Machine Learning, sect. 4, TR, 2004.
% [2] Jylanki, Vanhatalo & Vehtari, Robust GPR with a Student's-t
%     Likelihood, JMLR, 2011.
% [3] Bui, Yan & Turner, A Unifying Framework for Sparse GP Approximation
%     using Power EP, 2016, https://arxiv.org/abs/1605.07066
% [4] Snelson & Ghahramani, Sparse GPs using pseudo-inputs, NIPS, 2006.
% [5] Titsias, Var. Learning of Inducing Variables in Sparse GPs, AISTATS, 2009
% [6] Wilson & Nickisch, Kernel Interp. for Scalable Structured GPs, ICML, 2015
% [7] Han, Malioutov & Shin,  Large-scale Log-det Computation through Stochastic
%     Chebyshev Expansions, ICML, 2015.
% [8] Grigorievskiy, Lawrence, Sarkka, Parallelizable sparse inverse formulation
%     Gaussian processes (SpInGP),https://arxiv.org/pdf/1610.08035, 2016.
% [9] Sarkka, Solin, Hartikainen, Spatiotemporal learning via infinite-dim
%     Bayesian filtering and smoothing. IEEE Signal Processing Magazine, 2013.
%
% Copyright (c) by Carl Edward Rasmussen, Kun Dong, Insu Han and
%                                                    Hannes Nickisch 2018-06-03.
%
% See also cov/apxSparse.m, cov/apxGrid.m, gp.m.

if nargin<4, opt = []; end                           % make sure variable exists
if isnumeric(cov), c1 = 'numeric'; else c1 = cov{1}; end         % detect matrix
if isa(c1, 'function_handle'), c1 = func2str(c1); end         % turn into string
spars = strcmp(c1,'apxSparse') || strcmp(c1,'covFITC');
grid  = strcmp(c1,'apxGrid')   || strcmp(c1,'covGrid');
state = strcmp(c1,'apxState');
exact = ~grid && ~spars && ~state;

if exact                   % A) Exact computations using dense matrix operations
  if strcmp(c1,'numeric'), K = cov; dK = [];           % catch dense matrix case
  else
    [K,dK] = feval(cov{:},hyp.cov,x);     % covariance matrix and dir derivative
  end
  K = struct('mvm',@(x)mvmK_exact(K,x), 'P',@(x)x, 'Pt',@(x)x,... % mvm and proj 
             'fun',@(W)ldB2_exact(W,K,dK));

elseif spars                                          % B) Sparse approximations
  if isfield(opt,'s'), s = opt.s; else s = 1.0; end            % default is FITC
  xud = isfield(hyp,'xu');      % flag deciding whether to compute hyp.xu derivs
  if xud, cov{3} = hyp.xu; end                 % hyp.xu provided, replace cov{3}
  xu = cov{3}; nu = size(xu,1);                        % extract inducing points
  [Kuu,   dKuu]   = feval(cov{2}{:}, hyp.cov, xu);     % get the building blocks
  [Ku,    dKu]    = feval(cov{2}{:}, hyp.cov, xu, x);
  [diagK, ddiagK] = feval(cov{2}{:}, hyp.cov, x, 'diag');
  snud = isfield(hyp,'snu');   % flag deciding whether to compute hyp.snu derivs
  if snud, snu2 = exp(2*hyp.snu);                     % hyp.snu already provided
  else snu2 = 1e-6*(trace(Kuu)/nu);            % stabilise by 0.1% of signal std
  end
  Luu  = chol(Kuu+diag(snu2.*ones(nu,1)));         % Kuu + diag(snu2) = Luu'*Luu
  V  = Luu'\Ku;                                   % V = inv(Luu')*Ku => V'*V = Q
  g = max(diagK-sum(V.*V,1)',0);                         % g = diag(K) - diag(Q)
  K.mvm = @(x) V'*(V*x) + bsxfun(@times,s*g,x);   % efficient matrix-vector mult
  K.P = @(x) Luu\(V*x); K.Pt = @(x) V'*(Luu'\x);         % projection operations
  K.fun = @(W) ldB2_sparse(W,V,g,Luu,dKuu,dKu,ddiagK,s,xud,snud,snu2);

elseif grid                                            % C)  Grid approximations
  n = size(x,1);
  if isfield(opt,'cg_tol'), cgtol = opt.cg_tol;       % stop conjugate gradients
  else cgtol = 1e-6; end                                        % same as in pcg
  if isfield(opt,'cg_maxit'), cgmit = opt.cg_maxit;    % number of cg iterations
  else cgmit = min(n,20); end                                   % same as in pcg
  if isfield(opt,'deg'), deg = opt.deg; else deg = 3; end     % interpol. degree
  if isfield(opt,'stat'), stat = opt.stat; else stat = false; end    % show stat
  cgpar = {cgtol,cgmit}; xg = cov{3}; p = numel(xg);  % conjugate gradient, grid
  if isfield(opt,'ldB2_method'), meth=opt.ldB2_method; else meth='scale'; end
  m = 10; d = 15; mit = 50; sd = [];                    % set default parameters
  if strncmpi(meth,'cheby',5) || strncmpi(meth,'lancz',5)
    if isfield(opt,'ldB2_hutch'),        m =   opt.ldB2_hutch;     end
    if isfield(opt,'ldB2_cheby_degree'), d =   opt.ldB2_cheby_degree;    end
    if isfield(opt,'ldB2_maxit'),        mit = opt.ldB2_maxit;     end
    if isfield(opt,'ldB2_seed'),         sd =  opt.ldB2_seed;      end
    if stat
      fprintf('Stochastic (%s) logdet estimation (%d,%d,[d=%d])\n',meth,m,mit,d)
    end
  else
    if stat, fprintf('Scaled eigval logdet estimation\n'); end
  end
  ldpar = {meth,m,d,mit,[],sd};                              % logdet parameters
  [Kg,Mx] = feval(cov{:},hyp.cov,x,[],deg);  % grid cov structure, interp matrix
  if stat    % show some information about the nature of the p Kronecker factors
    fprintf(apxGrid('info',Kg,Mx,xg,deg));
  end
  K.mvm = @(x) Mx*Kg.mvm(Mx'*x);                    % mvm with covariance matrix
  K.P = @(x)x; K.Pt = @(x)x;                             % projection operations
  K.fun = @(W) ldB2_grid(W,K,Kg,xg,Mx,cgpar,ldpar); K.Mx = Mx;

elseif state
  if isfield(opt,'deriv'), deriv = opt.deriv; else deriv = true; end   % default
  if isfield(opt,'ngrid'), ng = max(opt.ngrid,6); else  ng = []; end   % default
  if isfield(opt,'balance'), bal = opt.balance; else  bal = false; end % default
  if isfield(opt,'qridge'),  qr = opt.qridge;   else  qr = 1e-7; end   % default
  if isfield(opt,'rridge'),  rr = opt.rridge;   else  rr = 0; end      % default
  if isfield(opt,'alpha_alg'), aalg = opt.alpha_alg;  else aalg = 'spingp'; end

  global t0, t0 = min(x);     % ugly but the only way to get this into covLINiso
  [F,L,Qc,H,Pinf,stat,dF,dQc,dPinf] = apxState(cov{2},hyp.cov);    % state space

  if bal            % balance state space model for improved numerical stability
    [T,F] = balance(F); L = T\L; H = H*T;                     % balance F,L,Qc,H
    LL = T\chol(Pinf,'lower'); Pinf = LL*LL';                     % balance Pinf
    for j=1:size(dF,3)                                    % balance dF and dPinf
      dF(:,:,j) = T\dF(:,:,j)*T; dPinf(:,:,j) = T\dPinf(:,:,j)/T;
    end
  end

  par = {'trans',F,L,Qc,H,Pinf,stat,x};
  if deriv
    [As,Qs,t,tid,dAs,dQs] = apxState(par{:},ng,dF,dQc,dPinf);
  else
    [As,Qs,t,tid] = apxState(par{:},ng); dAs = []; dQs = [];
  end
  d = numel(H); n = size(x,1); nh = size(dF,3);

  % construct G,Q,B matrices K = G*A*Q*A'*G', A = inv(B)
  G = [sparse(n,d), kron(sparse(tid(2:n+1),1:n,1,n,n),H)];      % contains index
  for i=1:n+1, Qs(:,:,i) = Qs(:,:,i) + qr*eye(d); end             % regularise Q
  Q = blktridiag([],Qs);
  B = blktridiag(-As(:,:,2:n+1),repmat(eye(d),[1,1,n+1]));
  dQ = cell(nh,1); dB = cell(nh,1);
  if deriv
    for i=1:nh
      dQi = reshape(dQs(:,:,i,:),d,d,n+1);
      dQ{i} = blktridiag([],dQi);
      dAi = reshape(dAs(:,:,i,2:n+1),d,d,n);
      dB{i} = blktridiag(-dAi,repmat(zeros(d),[1,1,n+1]));
    end
  end
  iQs = zeros(size(Qs)); for i=1:n+1, iQs(:,:,i) = inv(Qs(:,:,i)); end  % inv(Q)
  iQ = blktridiag([],iQs);
  K = struct('mvm',@(x)mvmK_state(x,G,Q,B), ... % mvm
             'P',@(x)x, 'Pt',@(x)x,...          % proj 
             'fun',@(W)ldB2_state(W,G,B,Q,iQ,As,Qs,H,Pinf,t,tid,rr,aalg,...
                                                         dB,dQ, dAs,dQs,dPinf));
end

%% A) Exact computations using dense matrix operations =========================
function [ldB2,solveKiW,dW,dldB2,L,triB] = ldB2_exact(W,K,dK)
  isWneg = any(W<0); n = numel(W);
  if isWneg                  % switch between Cholesky and LU decomposition mode
    A = eye(n) + bsxfun(@times,K,W');                     % asymmetric A = I+K*W
    [L,U,P] = lu(A); u = diag(U);         % compute LU decomposition, A = P'*L*U
    signU = prod(sign(u));                                           % sign of U
    detP = 1;               % compute sign (and det) of the permutation matrix P
    p = P*(1:n)';
    for i=1:n                                                     % swap entries
      if i~=p(i), detP = -detP; j = find(p==i); p([i,j]) = p([j,i]); end
    end
    if signU~=detP     % log becomes complex for negative values, encoded by inf
      ldB2 = Inf;
    else          % det(L) = 1 and U triangular => det(A) = det(P)*prod(diag(U))
      ldB2 = sum(log(abs(u)))/2;
    end                                            % compute inverse if required
    if nargout>1, Q = U\(L\P); solveKiW = @(r) bsxfun(@times,W,Q*r); end
    if nargout>4, L = -diag(W)*Q; end                              % overwrite L
  else                                                 % symmetric B = I+sW*K*sW
    sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K);             % Cholesky factor of B
    ldB2 = sum(log(diag(L)));                                    % log(det(B))/2
    solveKiW = @(r) bsxfun(@times,solve_chol(L,bsxfun(@times,r,sW)),sW);
    if nargout>2, Q = bsxfun(@times,1./sW,solve_chol(L,diag(sW))); end
  end
  if nargout>2
    dW = sum(Q.*K,2)/2;            % d log(det(B))/2 / d W = diag(inv(inv(K)+W))
    triB = trace(Q);                                      % triB = trace(inv(B))
    dldB2 = @(varargin) ldB2_deriv_exact(W,dK,Q, varargin{:});     % derivatives
  end

function dhyp = ldB2_deriv_exact(W,dK,Q, alpha,a,b)
  if nargin>3, R = alpha*alpha'; else R = 0; end
  if nargin>5, R = R + 2*a*b'; end
  dhyp.cov = dK( bsxfun(@times,Q,W) - R )/2;

function z = mvmK_exact(K,x)
  if size(x,2)==size(x,1) && max(max( abs(x-eye(size(x))) ))<eps      % x=eye(n)
    z = K;                             % avoid O(n^3) operation as it is trivial
  else
    z = K*x;
  end

%% B) Sparse approximations ====================================================
function [ldB2,solveKiW,dW,dldB2,L,triB] = ldB2_sparse(W,V,g,Luu,...
                                                dKuu,dKu,ddiagK,s,xud,snud,snu2)
  z = s*g.*W; t = 1/s*log(z+1); i = z<1e-4;  % s=0: t = g*W, s=1: t = log(g*W+1)
  t(i) = g(i).*W(i).*(1-z(i)/2+z(i).^2/3);         % 2nd order Taylor for tiny z
  dt = 1./(z+1); d = W.*dt;                               % evaluate derivatives
  nu = size(Luu,1); Vd = bsxfun(@times,V,d');
  Lu = chol(eye(nu) + V*Vd'); LuV = Lu'\V;               % Lu'*Lu=I+V*diag(d)*V'
  ldB2 = sum(log(diag(Lu))) + sum(t)/2;    % s=1 => t=log(g.*W+1), s=0 => t=g.*W
  md = @(r) bsxfun(@times,d,r); solveKiW = @(r) md(r) - md(LuV'*(LuV*md(r)));
  if nargout>2                % dW = d log(det(B))/2 / d W = diag(inv(inv(K)+W))
    dW = sum(LuV.*((LuV*Vd')*V),1)' + s*g.*d.*sum(LuV.*LuV,1)';
    dW = dt.*(g+sum(V.*V,1)'-dW)/2;                % add trace "correction" term
    dldB2 = @(varargin) ldB2_deriv_sparse(V,Luu,d,LuV,dKuu,dKu,ddiagK,s,...
                                                     xud,snud,snu2,varargin{:});
    if nargout>4
      L = solve_chol(Lu*Luu,eye(nu))-solve_chol(Luu,eye(nu));   % Sigma-inv(Kuu)
    end
    if nargout>5
      r = 1./(z+1); R = (eye(nu) + V*bsxfun(@times,V',r.*W))\V;
      triB = r'*(1-W.*sum(V.*R)'.*r);
    end
  end

function dhyp = ldB2_deriv_sparse(V,Luu,d,LuV,dKuu,dKu,ddiagK,s,...
                                                        xud,snud,snu2,alpha,a,b)
  % K + 1./W = V'*V + inv(D), D = diag(d)
  % Q = inv(K+inv(W)) = inv(V'*V + diag(1./d)) = diag(d) - LuVd'*LuVd;
  LuVd = bsxfun(@times,LuV,d'); diagQ = d - sum(LuVd.*LuVd,1)';
  F = Luu\V; Qu = bsxfun(@times,F,d') - (F*LuVd')*LuVd;
  if nargin>11, diagQ = diagQ-alpha.*alpha; Qu = Qu-(F*alpha)*alpha'; end
  Quu = Qu*F';
  if nargin>13
    diagQ = diagQ-2*a.*b; Qu = Qu-(F*a)*b'-(F*b)*a'; Quu = Quu-2*(F*a)*(F*b)';
  end
  diagQ = s*diagQ + (1-s)*d;                          % take care of s parameter
  Qu = Qu - bsxfun(@times,F,diagQ'); Quu = Quu - bsxfun(@times,F,diagQ')*F';
  if snud                                    % compute inducing noise derivative
    dhyp.snu = -diag(Quu).*snu2; if numel(snu2)==1, dhyp.snu=sum(dhyp.snu); end
  else nu = size(Quu,1); Quu = Quu + 1e-6*trace(Quu)/nu*eye(nu);   % fixed noise
  end
  if xud
    dhyp.cov = ddiagK(diagQ)/2; dhyp.xu = 0;
    [dc,dx] = dKu(Qu);   dhyp.cov = dhyp.cov + dc;   dhyp.xu = dhyp.xu + dx;
    [dc,dx] = dKuu(Quu); dhyp.cov = dhyp.cov - dc/2; dhyp.xu = dhyp.xu - dx/2;
  else
    dhyp.cov = ddiagK(diagQ)/2 + dKu(Qu) - dKuu(Quu)/2;
  end


%% C)  Grid approximations =====================================================
function [ldB2,solveKiW,dW,dldB2,L,triB] = ldB2_grid(W,K,Kg,xg,Mx,cgpar,ldpar)
  if all(W>=0)                                 % well-conditioned symmetric case
    sW = sqrt(W); msW = @(x) bsxfun(@times,sW,x);
    mvmB = @(x) msW(K.mvm(msW(x)))+x;
    solveKiW = @(r) msW(linsolve(msW(r),mvmB,cgpar{:}));
  else                 % less well-conditioned symmetric case if some negative W
    mvmKiW = @(x) K.mvm(x)+bsxfun(@times,1./W,x);
    solveKiW = @(r) linsolve(r,mvmKiW,cgpar{:});
  end                                                   % K*v = Mx*Kg.mvm(Mx'*v)
  dhyp.cov = [];                                                          % init
  if strncmpi(ldpar{1},'cheby',5)    % stochastic estim: logdet cheby/hutchinson
    dK = @(a,b) apxGrid('dirder',Kg,xg,Mx,a,b);
    if nargout<3            % save some computation depending on required output
      ldB2 = ldB2_cheby(W,K.mvm,dK, ldpar{2:end});
    else
      [ldB2,dhyp.cov,dW] = ldB2_cheby(W,K.mvm,dK, ldpar{2:end});
    end
  elseif strncmpi(ldpar{1},'lancz',5)% stochastic estim: logdet lancz/hutchinson
    % the sign of the maxit parameter encodes whether ARPACK is used or not
    mv = ver('matlab');     % get Matlab version, ARPACK interface not in Octave
    if numel(strfind(ldpar{1},'arpack'))==0 || numel(mv)==0
      ldpar{4} = -abs(ldpar{4});
    end
    dK = @(a,b) apxGrid('dirder',Kg,xg,Mx,a,b);
    if nargout<3            % save some computation depending on required output
      ldB2 = ldB2_lanczos(W,K.mvm,dK, ldpar{[2,4,6]});
    else    
      [ldB2,dhyp.cov,dW] = ldB2_lanczos(W,K.mvm,dK, ldpar{[2,4,6]});
    end
  else
    s = 3;                                    % Whittle embedding overlap factor
    [V,ee,e] = apxGrid('eigkron',Kg,xg,s);         % perform eigen-decomposition
    [ldB2,de,dW] = logdet_fiedler(e,W);     % Fiedler's upper bound, derivatives
    de = de.*double(e>0); % chain rule of max(e,0) in eigkron, Q = V*diag(de)*V'
    if nargout>3, dhyp.cov = ldB2_deriv_grid_fiedler(Kg,xg,V,ee,de,s); end
  end
  dldB2 = @(varargin) ldB2_deriv_grid(dhyp, Kg,xg,Mx, varargin{:});
  if ~isreal(ldB2), error('Too many negative W detected.'), end
  L = @(r) -K.P(solveKiW(K.Pt(r)));
  if nargout>5, error('triB not implemented'), end

function dhyp = ldB2_deriv_grid_fiedler(Kg,xg,V,ee,de,s)
  p = numel(Kg.kron);                              % number of Kronecker factors
  ng = [apxGrid('size',xg)',1];                                 % grid dimension
  dhyp = [];                              % dhyp(i) = trace( V*diag(de)*V'*dKi )
  for i=1:p
    z = reshape(de,ng); Vi = V{i};
    for j=1:p, if i~=j, z = apxGrid('tmul',ee{j}(:)',z,j); end, end
    if isnumeric(Vi)
      Q = bsxfun(@times,Vi,z(:)')*Vi';
      dhci = Kg.kron(i).dfactor(Q);
    else
      kii = Kg.kron(i).factor.kii;
      [junk,ni] = apxGrid('expand',xg{i}); di = numel(ni);
      xs = cell(di,1);                % generic (Strang) circular embedding grid
      for j=1:di, n2 = floor(ni(j)-1/2)+1; xs{j} = [1:n2,n2-2*ni(j)+2:0]'; end
      Fz = real(fftn(reshape(z,[ni(:)',1])));  % imaginary part of deriv is zero
      rep = 2*s*ones(1,di); if di==1, rep = [rep,1]; end           % replication
      Fzw = repmat(Fz,rep); % replicate i.e. perform transpose of circ operation
      [junk,xw] = apxGrid('circ',kii,ni,s);    % get Whittle circ embedding grid
      [junk,dkwi] = kii(apxGrid('expand',xw));             % evaluate derivative
      dhci = dkwi(Fzw(:));
    end
    dhyp = [dhyp; dhci];
  end

function dhyp = ldB2_deriv_grid(dhyp, Kg,xg,Mx, alpha,a,b)
  if nargin>4, dhyp.cov = dhyp.cov-apxGrid('dirder',Kg,xg,Mx,alpha,alpha)/2; end
  if nargin>6, dhyp.cov = dhyp.cov-apxGrid('dirder',Kg,xg,Mx,a,b);           end

function q = linsolve(p,mvm,varargin) % solve q = mvm(p) via conjugate gradients
  [q,flag,relres,iter] = conjgrad(mvm,p,varargin{:});                 % like pcg
  if ~flag,error('Not converged after %d iterations, r=%1.2e\n',iter,relres),end

% Solve x=A*b with symmetric A(n,n), b(n,m), x(n,m) using conjugate gradients.
% The method is along the lines of PCG but suited for matrix inputs b.
function [x,flag,relres,iter,r] = conjgrad(A,b,tol,maxit)
if nargin<3, tol = 1e-10; end
if nargin<4, maxit = min(size(b,1),20); end
x0 = zeros(size(b)); x = x0;
if isnumeric(A), r = b-A*x; else r = b-A(x); end, r2 = sum(r.*r,1); r2new = r2;
nb = sqrt(sum(b.*b,1)); flag = 0; iter = 1;
relres = sqrt(r2)./nb; todo = relres>=tol; if ~any(todo), flag = 1; return, end
on = ones(size(b,1),1); r = r(:,todo); d = r;
for iter = 2:maxit
  if isnumeric(A), z = A*d; else z = A(d); end
  a = r2(todo)./sum(d.*z,1);
  a = on*a;
  x(:,todo) = x(:,todo) + a.*d;
  r = r - a.*z;
  r2new(todo) = sum(r.*r,1);
  relres = sqrt(r2new)./nb; cnv = relres(todo)<tol; todo = relres>=tol;
  d = d(:,~cnv); r = r(:,~cnv);                           % get rid of converged
  if ~any(todo), flag = 1; return, end
  b = r2new./r2;                                               % Fletcher-Reeves
  d = r + (on*b(todo)).*d;
  r2 = r2new;
end

% Upper determinant bound on log |K*diag(W)+I| using Fiedler's 1971 inequality.
% K = kron( kron(...,K{2}), K{1} ), W = diag(w) both symmetric psd.
% The bound is exact for W = w*ones(N,1). Here, E = eig(K) are the
% eigenvalues of the matrix K.
% See also Prob.III.6.14 in Matrix Analysis, Bhatia 1997.
%
% Given nxn spd matrices C and D with ordered nx1 eigenvalues c, d 
% then det(C+D) <= prod(c+flipud(d))=exp(ub).
function [ub,dE,dW,ie,iw] = logdet_fiedler(E,W)
  [E,ie] = sort(E,'descend'); [W,iw] = sort(W,'descend');         % sort vectors
  N = numel(E); n = numel(W); k = n/N*E; % dimensions, approximate spectrum of K
  if n>N, k = [k;zeros(n-N,1)]; else k = k(1:n); end  % extend/shrink to match W
  kw1 = k.*W+1; ub = sum(log(kw1))/2;                     % evaluate upper bound
  dW = zeros(n,1); dW(iw) = k./(2*kw1); m = min(n,N);      % derivative w.r.t. W
  dE = zeros(N,1); dE(ie(1:m)) = W(1:m)./(N*2/n*kw1(1:m));  % deriative w.r.t. E

% Approximate l = log(det(B))/2, where B = I + sqrt(W)*K*sqrt(W) and compute
% the derivatives d l / d hyp w.r.t. covariance hyperparameters and
% d l / d W the gradient w.r.t. the precision matrix W.
%
% Large-scale Log-det Computation through Stochastic Chebyshev Expansions
% Insu Han, Dmitry Malioutov, Jinwoo Shin, ICML, 2015.
%
% Chebyshev polynomials T[0],..,T[d], where T[0](x)=1, T[1](x)=x, and
%                                           T[i+1](x)=2*x*T[i](x)-T[i-1](x).
% dT[0](x)=0, dT[1](x)=1, dT[i+1](x)=2*T[i](x)+2*x*dT[i](x)-dT[i-1](x)
%
% W      is the (diagonal) precision matrix represented by an nx1 vector
% K      is a function so that K(z) gives the mvm K*z
% dK     is a function so that dK(u,v) yields d trace(u'*K*v) / d hyp (dmvm = 2)
%                  or  so that dK(u)   yields d K*u / d hyp           (dmvm = 1)
%
% m      is the number of random probe vectors, default values is 10
% d      is the degree of the polynomial, default value is 15
%
% maxit  is the maximum number of iterations for the largest eigenvalue,
%                                                            default value is 50
% emax   is the maximum eigenvalue (in case it is known)
%
% seed   is the seed for generating the random probe vectors
% dmvm   indicates derivative type, default is 2
%
% Copyright (c) by Insu Han and Hannes Nickisch 2017-02-24.

function [ldB2,dhyp,dW] = ldB2_cheby(W,K,dK, m,d, maxit,emax, seed,dmvm)
  sW = sqrt(W); n = numel(W);
  B = @(x) x + bsxfun(@times,sW, K(bsxfun(@times,sW,x) ));%B=I+sqrt(W)*K*sqrt(W)
  if nargin<6, maxit = 50; end, if nargin<9, dmvm = 2; end        % set defaults
  if nargin<5, d = 15; end, if nargin<4, m = 10; end
  emaxconst = ~(nargin<7 || numel(emax)~=1);   % given as a constant => no deriv
  if ~emaxconst                                % evaluate upper eigenvalue bound
    if n==1, emax = B(1); else
      opt.maxit = maxit;   % number of iterations to estimate the largest eigval
      opt.issym = 1; opt.isreal = 1; % K is real symmetric and - of course - psd
      cstr = 'eigs(B,n,1,''lm'',opt)';              % compute largest eigenvalue
      if exist('evalc'), [txt,vmax,emax] = evalc(cstr);
      else [vmax,emax] = eval(cstr); end
    end
    if nargout>1                        % compute  d n*log(1+emax)/2 / d {hyp,W}
      r = sW.*vmax; dW = (K(r).*(r./W) + K(r).*(r./W))/2;       % d emax / d W
      if dmvm==1, dhyp = dK(r)'*r; else dhyp = dK(r,r); end     % d emax / d hyp
      dW   = n/(2*(1+emax)) * dW; dhyp = n/(2*(1+emax)) * dhyp;        % rescale
      % for the second term in ldB2 coming from the Chebyshev polynomial, we
      % ignore the dependency on the largest eigenvalue emax, which essentially
      % renders the derivative slightly incorrect
    end
  else
    dhyp = 0; dW = zeros(n,1);                      % init derivatives with zero
  end
  d = round(abs(d)); if d==0, error('We require d>0.'), end           % emin = 1
  a = 1+emax; del = 1-emax/a; ldB2 = n*log(a)/2;   % scale eig(B) to [del,1-del]
  if emax>1e5                                                   % del=1/(1+emax)
    fprintf('B has large condition number %1.2e\n',emax)
    fprintf('log(det(B))/2 is most likely overestimated\n')
  end
  s = 2/(emax-1);                    % compute scaling factor, s = 2/a/(1-2*del)
  A = @(x) s*B(x) - s*a/2*x;         % apply scaling transform: [1,emax]->[-1,1]

  xk = cos(pi*((0:d)'+0.5)/(d+1));                              % zeros of Tn(x)
  fk = log(((1-2*del)/2).*xk+0.5);        % target function, [-1,1]->[del,1-del]
  Tk = ones(d+1,d+1); Tk(:,2) = xk;                             % init recursion
  for i=2:d, Tk(:,i+1) = 2*xk.*Tk(:,i) - Tk(:,i-1); end   % evaluate polynomials
  c = 2/(d+1)*(Tk'*fk); c(1) = c(1)/2;          % compute Chebyshev coefficients
  if nargin>7 && numel(seed)>0, randn('seed',seed), end               % use seed
  V = sign(randn(n,m));                         % bulk draw Rademacher variables
  p1 = [1; zeros(d,1)]; p2 = [0;1;zeros(d-1,1)]; % Chebyshev->usual coefficients
  p = c(1)*p1 + c(2)*p2;
  for i=2:d, p3 = [0;2*p2(1:d)]-p1; p = p + c(i+1)*p3; p1 = p2; p2 = p3; end
  if nargout<2  % no derivs requested; use bulk MVMs with A, one for each j=1..d
%     U = p(1)*V; AjV = V; for j=1:d, AjV = A(AjV); U = U + p(j+1)*AjV; end
    U1 = V; U2 = A(V); U = c(1)*U1 + c(2)*U2;          % numerically more robust
    for i=2:d, U3 = 2*A(U2) - U1; U = U+c(i+1)*U3; U1 = U2; U2 = U3; end
    ldB2 = ldB2 + sum(sum(V.*U))/(2*m);    % sum_{j=0..d} p(j+1)*trace(V'*A^j*U)
  elseif dmvm==1                                     % use MVM-based derivatives
    dA = @(x) s*bsxfun(@times,sW, dK(bsxfun(@times,sW,x) ));    % derivative MVM
    if nargout>2
      if norm(diff(W))>1e-12, error('Only isotropic dW supported.'), end
      dA = @(x) [dA(x), s*K(x)];        % augment to get deriv w.r.t. W=w*eye(n)
    end
    U = p(1)*V; AjV = V; dU = 0;
    for j=1:d
      if j==1, dAjV = dA(AjV); else dAjV = dA(AjV) + A(dAjV); end
      dU = dU + p(j+1)*dAjV; AjV = A(AjV); U = U + p(j+1)*AjV;
    end
    nh = size(dAjV,2)/size(AjV,2); dhyp = zeros(nh,1);
    for j=1:nh, dhyp(j) = sum(sum(V.*dU(:,(j-1)*m+(1:m))))/(2*m); end
    ldB2 = ldB2 + sum(sum(V.*U))/(2*m);    % sum_{j=0..d} p(j+1)*trace(V'*B^j*U)
    if nargout>2, dW = dW + dhyp(nh)/n; nh = nh-1; dhyp = dhyp(1:nh); end
  else                                   % deal with one random vector at a time
    Av = zeros(n,d+1);                            % all powers Av(:,j) = (A^j)*v
    for i=1:m
      v = V(:,i); Av(:,1) = v;
      for j=1:d, Av(:,j+1) = A(Av(:,j)); end
      ldB2 = ldB2 + (v'*Av*p)/(2*m);
      sWAv = bsxfun(@times,sW,Av);
      for j=1:d                   % p(1)*I + p(2)*A + p(3)*A^2 + .. + p(d+1)*A^d
        akj = s*p(j+1)/m;
        % equivalent to: k = 1:j; dhyp = dhyp + akj*dK(sWAv(:,j-k+1),sWAv(:,k));
        % exploiting symmetry dK(a,b)=dK(b,a) reduces the computation to half
        k = 1:ceil((j-1)/2); dhyp = dhyp + akj  *dK(sWAv(:,j-k+1),sWAv(:,k));
        if mod(j,2)
          k = ceil(j/2); dhyp = dhyp + akj/2*dK(sWAv(:,j-k+1),sWAv(:,k));
        end
      end
      if nargout>2, u = bsxfun(@times,1./sW,Av); w = K(sWAv);       % precompute
        for j=1:d
          dWj = sum( u(:,j:-1:1).*w(:,1:j) + w(:,j:-1:1).*u(:,1:j), 2 );
          dW = dW + s*p(j+1)/(4*m) * dWj;
        end
      end
    end
  end

% Approximate l = log(det(B))/2, where B = I + sqrt(W)*K*sqrt(W) and
% the derivatives d l / d hyp w.r.t. covariance hyperparameters and
% d l / d W the gradient w.r.t. the (effective) precision matrix W.
%
% W      is the (diagonal) precision matrix represented by an nx1 vector
% K      is a function so that K(z) gives the mvm K*z
% dK     is a function so that dK(u,v) yields d trace(u'*K*v) / d hyp (dmvm = 2)
%                  or  so that dK(u)   yields d K*u / d hyp           (dmvm = 1)
%
% m      is the number of random probe vectors, default values is 10
% maxit  is the number of Lanczos vectors, default value is 15
%
% seed   is the seed for generating the random probe vectors
% dmvm   indicates derivative type, default is 2
%
% Copyright (c) by Kun Dong and Hannes Nickisch 2017-05-05.

function [ld,dhyp,dW] = ldB2_lanczos(W,K,dK, m,maxit, seed,dmvm)
  n = numel(W); sW = sqrt(W); if nargin<6, seed = 42; end   % massage parameters
  if nargin<5, maxit = 15; end, if nargin<4, m = 10; end    % set default values
  arpack = maxit>0; maxit = abs(maxit);       % sign encodes ARPACK versus plain
  if size(m,2)>1, V = m(1:n,:); m = size(m,2); else V = sign(randn(n,m)); end
  if nargin>5 && numel(seed)>0, randn('seed',seed), end               % use seed
  if nargin<7, dmvm = 2; end, ld = 0; dhyp = 0; dW = 0;           % set defaults
  sWm = @(v) bsxfun(@times,sW,v);                            % MVM with diag(sW)
  B = @(v) v + sWm(K(sWm(v)));
  for j = 1:m
    if arpack
      [Q,T] = infGrid('lanczos_arpack',B,V(:,j),maxit);    % Lanczos with ARPACK
    else
      [Q,T] = lanczos_full(B,V(:,j),maxit);              % perform plain Lanczos
    end
    [Y,f] = eig(T); f = diag(f);
    ld = ld + n/(2*m) * (Y(1,:).^2*log(f));        % evaluate Stieltjes integral
    if nargout > 1                                        % go after derivatives
      v0 = sW.*V(:,j)/norm(V(:,j));                        % avoid call to normc
      w0 = (Q*(T\[1; zeros(numel(f)-1,1)])).*sW;   
      if dmvm==1, dhypj = dK(v0)'*w0;  else dhypj = dK(v0,w0); end
      dhyp = dhyp + n/(2*m) * dhypj;
      if nargout>2, dW = dW + n/(2*m)*(K(v0).*(w0./W) + K(w0).*(v0./W))/2; end
    end
  end

function [Q,T] = lanczos_full(B,v,d)       % perform Lanczos with at most d MVMs
  alpha = zeros(d,1); beta = zeros(d,1);                                 % for T
  v = v/norm(v,'fro'); n = numel(v);
  u = zeros(n,1); k = 1; Q = zeros(n,d);                   % mem for alg and res
  for k=1:d                                            % do at most d iterations
    [u,v,alpha(k),beta(k)] = lanczos_step(u,v,B,Q(:,1:k-1));
    Q(:,k) = u;                          % store Lanczos vectors for later usage
    if abs(beta(k))<1e-10, break, end                     % diagnose convergence
  end
  alpha = alpha(1:k); beta = beta(2:k); Q = Q(:,1:k);                 % truncate
  T = diag(alpha) + diag(beta,1 )+ diag(beta,-1);

function [u,v,a,b] = lanczos_step(u,v,B,Q)       % Lanczos step, $9.2.1 of Golub
  b = norm(v,'fro'); t = u; u = v/b;
  u = u - Q*(Q'*u); u = u/norm(u);                        % perform Gram-Schmidt
  r = B(u)-b*t; a = real(u'*r); v = r-a*u;


%% D)  State space approximations ==============================================
function [ldB2,solveKiW,dW,dldB2,L,triB] = ldB2_state(W,G,B,Q,iQ,As,Qs,...
                                      H,Pinf,t,tid,rr,aalg,dB,dQ, dAs,dQs,dPinf)
  n = numel(W); s = apxState('kalman',As,Qs, H,Pinf, zeros(n,1),W, t,tid);
  ldB2 = sum(log(s))/2; w = W;          % ldB2 = log(det(I+sqrt(W)*K*sqrt(W)))/2
  % inv(K+inv(diag(W))) = diag(W) - diag(W)*G*inv(R)*G'*diag(W) where R is BTD
  W = sparse(1:n,1:n,W,n,n); R = B'*iQ*B + G'*W*G;    % R = inv(Kt)+G'*diag(W)*G
  R = (R+R')/2;                                               % enforce symmetry
  if rr>0, R = R + rr*speye(size(R)); end                           % regularize

  if strcmp(aalg,'spingp')  % 2 options: either through sparse matrix algebra ..
    solveKiW = @(r) W*r - W*(G*(R\(G'*(W*r))));
  else                                    % .. or via Kalman filtering+smoothing
    solveKiW = @(r) solveKiW_state(As,Qs, H,Pinf, w.*r,w, t,tid);
  end                                           % dW = diag(G*inv(R)*G')/2 = ..
  if nargout<=2, return, end             % .. diag(inv(I+sqrt(W)*K*sqrt(W))*K)/2
  if exist('sparseinv_mex','file')
    iR = sparseinv(R);                                    % scalable computation
  else
    iR = R\speye(size(R));                    % this computation is NOT scalable
  end
  dW = full(diag(G*iR*G'))/2;
  nh = numel(dB); z = zeros(n,1);
  if nh>0 && numel(dB{1})>0
    [s,ga,nu,tau,ms,Ps,ds] = apxState('kalman',As,Qs, H,Pinf, z,w,t,tid,[], ...
                                                                 dAs,dQs,dPinf);
    tr = ds'*(1./s)/2;% tr(i)=tr(dKi*Z),Z=W-W*(G*(R\(G'*W)))=inv(K+inv(diag(W)))
  else
    tr = zeros(nh,1);
  end
  dldB2 = @(varargin) ldB2_deriv_state(G,B,dB,Q,dQ,tr, varargin{:});
  L = @(r) -solveKiW(r);
  triB = n - 2*sum(W*dW);               % triB = trace(inv(I+sqrt(W)*K*sqrt(W)))

function al = solveKiW_state(As,Qs, H,Pinf, r,W, t,tid)
  [z,ga, tnu,ttau, ms,Ps] = apxState('kalman',As,Qs, H,Pinf, r,W, t,tid);
  [ms,Ps,rho] = apxState('rts',As,Qs,ms,Ps,H);
  rho(tid) = rho; al = ga - W.*rho(1:end-1);

% dhyp(i) = tr( (Q-R)'*dKi )/2, where dKi = d K / d hyp(i)
% R = alpha*alpha' + 2*a*b', Q = inv(K+inv(W))
function dhyp = ldB2_deriv_state(G,B,dB,Q,dQ,tr, alpha,a,b)
  nh = numel(dB); dhyp.cov = zeros(nh,1);  % number of hyperparameters, allocate
  if nh==0 || numel(dB{1})==0, return, end    % no matrices available, stop here
  for i=1:nh %                           iterate over covariance hyperparameters
    dKi = @(x) dmvmK_state(x,G,Q,B,dQ{i},dB{i});       % MVM with d K / d hyp(i)
    dhyp.cov(i) = tr(i);
    if nargin>6, dhyp.cov(i) = dhyp.cov(i) - alpha'*dKi(alpha)/2; end
    if nargin>8, dhyp.cov(i) = dhyp.cov(i) - a'*dKi(b); end
  end

function z = dmvmK_state(x,G,Q,B,dQ,dB)      % matrix-vector multiplication O(n)
  a = B'\(G'*x); z = G*(B\(dQ*a))-G*(B\(dB*(B\(Q*a))))-G*(B\(Q*(B'\(dB'*a))));

function z = mvmK_state(x,G,Q,B)             % matrix-vector multiplication O(n)
  z = G*(B\(Q*(B'\(G'*x))));

% Construct and solve the block tridiagonal system A*b = x with A of shape
% (q*n,q*n) and blocksize q.
% For q=1, we have A = diag(l(:),-1) + diag(d(:)) + diag(u(:).
%
%  x and b have shape (q*n,m)
%  l and u have shape (n-1,1) for q=1 or (q,q,n-1) or (0,0)
%  d       has  shape (n,1)   for q=1 or (q,q,n)
%
%  [ d(1)  u(1)                                  ] [  b(1)  ]   [  x(1)  ]
%  [ l(1)  d(2)  u(2)                            ] [  b(2)  ]   [  x(2)  ]
%  [       l(2)  d(3)  u(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        l(n-2) d(n-1) u(n-1) ] [ b(n-1) ]   [ x(n-1) ]
%  [                               l(n-1)  u(n)  ] [  b(n)  ]   [  x(n)  ]
function [A,b] = blktridiag(l,d,u, x)
  if ndims(l)==2 && numel(l)>0, l = reshape(l,1,1,[]); end
  if ndims(d)==2, d = reshape(d,1,1,[]); end, q = size(d,1); n = size(d,3);
  if nargin<3, u = []; end                                 % default is u(:) = 0
  if ndims(u)==2 && numel(u)>0, u = reshape(u,1,1,[]); end

  o = reshape(repmat(1:q*n,q,1),q,q,n); vec = @(x) x(:);
  ii = vec(permute(o,[2,1,3])); jj = vec(o);      % block diagonal index vectors
  i = ii; j = jj; s = vec(d);
  if numel(l)>0, i=[i;ii(q*q+1:end)  ]; j=[j;jj(q*q+1:end)-q]; s=[s;vec(l)]; end
  if numel(u)>0, i=[i;ii(q*q+1:end)-q]; j=[j;jj(q*q+1:end)  ]; s=[s;vec(u)]; end
  A = sparse(i,j,s,q*n,q*n);
  if nargout>1, b = A\x; end                     % only solve system if required