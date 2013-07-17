function [post nlZ dnlZ] = infLaplace(hyp, mean, cov, lik, x, y)

% Laplace approximation to the posterior Gaussian process.
% The function takes a specified covariance function (see covFunction.m) and
% likelihood function (see likFunction.m), and is designed to be used with
% gp.m. See also infFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2012-11-09.
%
% See also INFMETHODS.M.

persistent last_alpha                                   % copy of the last alpha
if any(isnan(last_alpha)), last_alpha = zeros(size(last_alpha)); end   % prevent
tol = 1e-6;                   % tolerance for when to stop the Newton iterations
smax = 2; Nline = 10; thr = 1e-4;                       % line search parameters

maxit = 20;                                    % max number of Newton steps in f

inf = 'infLaplace';
n = size(x,1);
K = feval(cov{:},  hyp.cov,  x);                % evaluate the covariance matrix
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector

Psi_old = Inf;    % make sure while loop starts by the largest old objective val
if any(size(last_alpha)~=[n,1])     % find a good starting point for alpha and f
  alpha = zeros(n,1); f = K*alpha+m;       % start at mean if sizes do not match                  
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W = -d2lp; Psi_new=-sum(lp);
else
  alpha = last_alpha; f = K*alpha+m;                              % try last one
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W = -d2lp;
  Psi_new = alpha'*(f-m)/2 - sum(lp);                 % objective for last alpha
  Psi_def = -sum(feval(lik{:},hyp.lik,y,m,[],inf));% objective for def init f==m
  if Psi_def < Psi_new                         % if default is better, we use it
    alpha = zeros(n,1); f = K*alpha+m;
    [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp; Psi_new=-sum(lp);
  end
end
isWneg = any(W<0);       % flag indicating whether we found negative values of W
it = 0;                            % this happens for the Student's t likelihood
while Psi_old - Psi_new > tol && it<maxit                         % begin Newton
  Psi_old = Psi_new; it = it+1;
  if isWneg       % stabilise the Newton direction in case W has negative values
    W = max(W,0);      % stabilise the Hessian to guarantee postive definiteness
    tol = 1e-10;           % increase accuracy to also get the derivatives right
    % In Vanhatalo et. al., GPR with Student's t likelihood, NIPS 2009, they use
    % a more conservative strategy then we do being equivalent to 2 lines below.
    % nu  = exp(hyp.lik(1));                  % degree of freedom hyperparameter
    % W  = W + 2/(nu+1)*dlp.^2;               % add ridge according to Vanhatalo
  end
  sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K);              % L'*L=B=eye(n)+sW*K*sW
  b = W.*(f-m) + dlp;
  dalpha = b - sW.*solve_chol(L,sW.*(K*b)) - alpha;   % Newton dir + line search
  [s,Psi_new,Nfun,alpha,f,dlp,W] = brentmin(0,smax,Nline,thr, ...
                                    @Psi_line,4,dalpha,alpha,hyp,K,m,lik,y,inf);
  isWneg = any(W<0);
end                                                    % end Newton's iterations

last_alpha = alpha;                                     % remember for next call
[lp,dlp,d2lp,d3lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp; isWneg=any(W<0);
post.alpha = alpha;                            % return the posterior parameters
post.sW = sqrt(abs(W)).*sign(W);             % preserve sign in case of negative

if isWneg                    % switch between Cholesky and LU decomposition mode
  % For post.L = -inv(K+diag(1./W)), we us the non-default parametrisation.
  [ldA, iA, post.L] = logdetA(K,W);     % A=eye(n)+K*W is as safe as symmetric B
  nlZ = alpha'*(f-m)/2 - sum(lp) + ldA/2;
else
  sW = post.sW; post.L = chol(eye(n)+sW*sW'.*K);                     % recompute
  nlZ = alpha'*(f-m)/2 + sum(log(diag(post.L))-lp);     % ..(f-m)/2 -lp +ln|B|/2
end

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  if isWneg                  % switch between Cholesky and LU decomposition mode
    Z = -post.L;                                                 % inv(K+inv(W))
    g = sum(iA.*K,2)/2; % deriv. of ln|B| wrt W; g = diag(inv(inv(K)+diag(W)))/2
  else
    Z = repmat(sW,1,n).*solve_chol(post.L,diag(sW)); %sW*inv(B)*sW=inv(K+inv(W))
    C = post.L'\(repmat(sW,1,n).*K);                     % deriv. of ln|B| wrt W
    g = (diag(K)-sum(C.^2,1)')/2;                    % g = diag(inv(inv(K)+W))/2
  end
  dfhat = g.*d3lp;  % deriv. of nlZ wrt. fhat: dfhat=diag(inv(inv(K)+W)).*d3lp/2
  for i=1:length(hyp.cov)                                    % covariance hypers
    dK = feval(cov{:}, hyp.cov, x, [], i);
    dnlZ.cov(i) = sum(sum(Z.*dK))/2 - alpha'*dK*alpha/2;         % explicit part
    b = dK*dlp;                            % b-K*(Z*b) = inv(eye(n)+K*diag(W))*b
    dnlZ.cov(i) = dnlZ.cov(i) - dfhat'*( b-K*(Z*b) );            % implicit part
  end
  for i=1:length(hyp.lik)                                    % likelihood hypers
    [lp_dhyp,dlp_dhyp,d2lp_dhyp] = feval(lik{:},hyp.lik,y,f,[],inf,i);
    dnlZ.lik(i) = -g'*d2lp_dhyp - sum(lp_dhyp);                  % explicit part
    b = K*dlp_dhyp;                        % b-K*(Z*b) = inv(eye(n)+K*diag(W))*b
    dnlZ.lik(i) = dnlZ.lik(i) - dfhat'*( b-K*(Z*b) );            % implicit part
  end
  for i=1:length(hyp.mean)                                         % mean hypers
    dm = feval(mean{:}, hyp.mean, x, i);
    dnlZ.mean(i) = -alpha'*dm;                                   % explicit part
    dnlZ.mean(i) = dnlZ.mean(i) - dfhat'*(dm-K*(Z*dm));          % implicit part
  end
end

% criterion Psi at alpha + s*dalpha for line search
function [Psi,alpha,f,dlp,W] = Psi_line(s,dalpha,alpha,hyp,K,m,lik,y,inf)
  alpha = alpha + s*dalpha;
  f = K*alpha+m;
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W = -d2lp;
  Psi = alpha'*(f-m)/2 - sum(lp);

% Compute the log determinant ldA and the inverse iA of a square nxn matrix
% A = eye(n) + K*diag(w) from its LU decomposition; for negative definite A, we 
% return ldA = Inf. We also return mwiA = -diag(w)*inv(A).
function [ldA,iA,mwiA] = logdetA(K,w)
  [m,n] = size(K); if m~=n, error('K has to be nxn'), end
  A = eye(n)+K.*repmat(w',n,1);
  [L,U,P] = lu(A); u = diag(U);           % compute LU decomposition, A = P'*L*U
  signU = prod(sign(u));                                             % sign of U
  detP = 1;                 % compute sign (and det) of the permutation matrix P
  p = P*(1:n)';
  for i=1:n                                                       % swap entries
    if i~=p(i), detP = -detP; j = find(p==i); p([i,j]) = p([j,i]); end
  end
  if signU~=detP  % log becomes complex for negative values, encoded by infinity
    ldA = Inf;
  else            % det(L) = 1 and U triangular => det(A) = det(P)*prod(diag(U))
    ldA = sum(log(abs(u)));
  end 
  if nargout>1, iA = inv(U)*inv(L)*P; end          % return the inverse, as well
  if nargout>2, mwiA = -repmat(w,1,n).*iA; end
