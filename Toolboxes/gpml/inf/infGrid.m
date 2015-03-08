function [post nlZ dnlZ] = infGrid(hyp, mean, cov, lik, x, y, opt)

% Inference for a GP with Gaussian likelihood and covGrid covariance.
% The (Kronecker) covariance matrix used is given by:
%   K = kron( kron(...,K{2}), K{1} ) = K_p x .. x K_2 x K_1.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters.
% The result is exact for complete grids, otherwise results are approximate.
% See also "help infMethods".
%
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with covFITC and likGauss.
%
% Copyright (c) by Hannes Nickisch and Andrew Wilson, 2014-12-03.
%
% See also INFMETHODS.M, COVGRID.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Inference with infGrid only possible with Gaussian likelihood.');
end
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if ~strcmp(cov1,'covGrid'); error('Only covGrid supported.'), end    % check cov

xg = cov{3}; p = numel(xg);                                    % underlying grid
N = 1; D = 0; for i=1:p, N = N*size(xg{i},1); D = D+size(xg{i},2); end    % dims
[K,M,xe] = feval(cov{:}, hyp.cov, x);     % evaluate covariance mat constituents
xe = M*xe; n = size(xe,1);
m = feval(mean{:}, hyp.mean, xe);                         % evaluate mean vector

if nargin<=6, opt = []; end                        % make opt variable available
if isfield(opt,'cg_tol'),   cgtol = opt.cg_tol;       % stop conjugate gradients
else cgtol = 1e-5; end
if isfield(opt,'cg_maxit'), cgmit = opt.cg_maxit;      % number of cg iterations
else cgmit = min(n,20); end
for j=1:numel(K)
  if iscell(K{j}) && strcmp(K{j}{1},'toep'), K{j} = toeplitz(K{j}{2}); end
end

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
V = cell(size(K)); E = cell(size(K));                 % eigenvalue decomposition
for j=1:numel(K)
  if iscell(K{j}) && strcmp(K{j}{1},'toep')
       [V{j},E{j}] = eigr(toeplitz(K{j}{2}));
  else [V{j},E{j}] = eigr(K{j});
  end
end
e = 1; for i=1:numel(E), e = kron(diag(E{i}),e); end       % eigenvalue diagonal
if numel(m)==N
  s = 1./(e+sn2); ord = 1:N;                  % V*diag(s)*V' = inv(K+sn2*eye(N))
else
  [eord,ord] = sort(e,'descend');              % sort long vector of eigenvalues
  s = 1./((n/N)*eord(1:n) + sn2);               % approx using top n eigenvalues
end

if numel(m)==N               % decide for Kronecker magic or conjugate gradients
  L = @(k) kronmvm(V,repmat(-s,1,size(k,2)).*kronmvm(V,k,1));     % mvm callback
else                                                % go for conjugate gradients
  mvm = @(t) mvmK(t,K,sn2,M);                             % single parameter mvm
  L = @(k) -solveMVM(k,mvm,cgtol,cgmit);
end
alpha = -L(y-m);

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                         % sqrt of noise precision
post.L = L;                           % function to compute inv(K+sn2*eye(N))*Ks

if nargout>1                               % do we want the marginal likelihood?
  lda = -sum(log(s));                                               % exact kron
  % exact: lda = 2*sum(log(diag(chol(M*kronmvm(K,eye(N))*M'+sn2*eye(n)))));
  nlZ = (y-m)'*alpha/2 + n*log(2*pi)/2 + lda/2;
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;     % allocate space for derivatives, define Q=inv(M*K*M'+sn2*I)
    Mtal = M'*alpha;                          % blow up alpha vector from n to N
    for i = 1:numel(hyp.cov)
      dK = feval(cov{:}, hyp.cov, x, [], i);
      for j=1:numel(dK)                                        % expand Toeplitz
        if iscell(dK{j}) && strcmp(dK{j}{1},'toep')
          dK{j} = toeplitz(dK{j}{2});
        end
      end
      P = cell(size(dK));               % auxiliary Kronecker matrix P = V'*dK*V
      for j = 1:numel(dK), P{j} = sum((V{j}'*dK{j}).*V{j}',2); end
      p = 1; for j=1:numel(P), p = kron(P{j},p); end               % p = diag(P)
      p = (n/N)*p(ord(1:n));        % approximate if incomplete grid observation
      dnlZ.cov(i) = (p'*s - Mtal'*kronmvm(dK,Mtal))/2;
    end
    dnlZ.lik = sn2*(sum(s) - alpha'*alpha);                  % sum(s) = trace(Q)
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, xe, i)'*alpha;
    end
  end
end

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

% solve q = mvm(p) via conjugate gradients
function q = solveMVM(p,mvm,varargin)
  [q,flag,relres,iter] = conjgrad(mvm,p,varargin{:});                 % like pcg
  if ~flag,error('Not converged after %d iterations, r=%1.2e\n',iter,relres),end

% mvm so that q = M*K*M'*p + sn2*p using the Kronecker representation
function q = mvmK(p,K,sn2,M)
  q = M*kronmvm(K,M'*p) + sn2*p;

% Perform a matrix vector multiplication b = A*x with a matrix A being a
% Kronecker product given by A = kron( kron(...,As{2}), As{1} ).
function b = kronmvm(As,x,transp)
if nargin>2 && ~isempty(transp) && transp   % transposition by transposing parts
  for i=1:numel(As), As{i} = As{i}'; end
end
m = zeros(numel(As),1); n = zeros(numel(As),1);                  % extract sizes
for i=1:numel(n)
  if iscell(As{i}) && strcmp(As{i}{1},'toep')
    m(i) = size(As{i}{2},1); n(i) = size(As{i}{2},1);
  else [m(i),n(i)] = size(As{i});
  end
end
d = size(x,2);
b = x;
for i=1:numel(n)
  a = reshape(b,[prod(m(1:i-1)), n(i), prod(n(i+1:end))*d]);    % prepare  input
  tmp = reshape(permute(a,[1,3,2]),[],n(i))*As{i}';
  b = permute(reshape(tmp,[size(a,1),size(a,3),m(i)]),[1,3,2]);
end
b = reshape(b,prod(m),d);                        % bring result in correct shape

% Real eigenvalues and eigenvectors up to the rank of a real symmetric matrix.
% Decompose A into V*D*V' with orthonormal matrix V and diagonal matrix D.
% Diagonal entries of D obave the rank r of the matrix A as returned by
% the call rank(A,tol) are zero.
function [V,D] = eigr(A,tol)
[V,D] = eig((A+A')/2); n = size(A,1);    % decomposition of strictly symmetric A
d = max(real(diag(D)),0); [d,ord] = sort(d,'descend');        % tidy up and sort
if nargin<2, tol = size(A,1)*eps(max(d)); end, r = sum(d>tol);     % get rank(A)
d(r+1:n) = 0; D = diag(d);                % set junk eigenvalues to strict zeros
V(:,1:r) = real(V(:,ord(1:r))); V(:,r+1:n) = null(V(:,1:r)'); % ortho completion
