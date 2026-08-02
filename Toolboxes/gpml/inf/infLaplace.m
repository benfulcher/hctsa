function [post nlZ dnlZ alpha] = infLaplace(hyp, mean, cov, lik, x, y, opt)

% Laplace approximation to the posterior Gaussian process.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-10-21.
%
% See also infMethods.m, gp.m.

persistent last_alpha                                   % copy of the last alpha
if any(isnan(last_alpha)), last_alpha = zeros(size(last_alpha)); end   % prevent
if nargin<=6, opt = []; end                        % make opt variable available

inf = 'infLaplace';
n = size(x,1);
if isstruct(cov), K = cov;                   % use provided covariance structure
else K = apx(hyp,cov,x,opt); end               % set up covariance approximation
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv
likfun = @(f) feval(lik{:},hyp.lik,y,f,[],inf);        % log likelihood function

if any(size(last_alpha)~=[n,1])     % find a good starting point for alpha and f
  alpha = zeros(n,1);                      % start at mean if sizes do not match
else
  alpha = last_alpha;                                             % try last one
  if Psi(alpha,m,K.mvm,likfun) > -sum(likfun(m)) % default f==m better => use it
    alpha = zeros(n,1);
  end
end

% switch between optimisation methods
alpha = irls(alpha, m,K,likfun, opt);                         % run optimisation

f = K.mvm(alpha)+m;                             % compute latent function values
last_alpha = alpha;                                     % remember for next call
[lp,dlp,d2lp,d3lp] = likfun(f); W = -d2lp;
[ldB2,solveKiW,dW,dhyp] = K.fun(W);        % obtain functionality depending on W
post.alpha = K.P(alpha);                       % return the posterior parameters
post.sW = sqrt(abs(W)).*sign(W);             % preserve sign in case of negative
post.L = @(r) -K.P(solveKiW(K.Pt(r)));

% diagnose optimality
err = @(x,y) norm(x-y)/max([norm(x),norm(y),1]);   % we need to have alpha = dlp
% dev = err(alpha,dlp);  if dev>1e-4, warning('Not at optimum %1.2e.',dev), end

nlZ = alpha'*(f-m)/2 - sum(lp) + ldB2;             % compute marginal likelihood
if nargout>2                                           % do we want derivatives?
  dfhat = dW.*d3lp; % deriv. of nlZ wrt. fhat: dfhat=diag(inv(inv(K)+W)).*d3lp/2
  dahat = dfhat - solveKiW(K.mvm(dfhat)); dnlZ = dhyp(alpha,dlp,dahat);
  dnlZ.lik = zeros(size(hyp.lik));                             % allocate memory
  for i=1:length(hyp.lik)                                    % likelihood hypers
    [lp_dhyp,dlp_dhyp,d2lp_dhyp] = feval(lik{:},hyp.lik,y,f,[],inf,i);
    dnlZ.lik(i) = -dW'*d2lp_dhyp - sum(lp_dhyp);                 % explicit part
    b = K.mvm(dlp_dhyp);                   % b-K*(Z*b) = inv(eye(n)+K*diag(W))*b
    dnlZ.lik(i) = dnlZ.lik(i) - dfhat'*( b-K.mvm(solveKiW(b)) );      % implicit
  end
  dnlZ.mean = -dm(alpha+dahat);                       % explicit + implicit part
end

% Evaluate criterion Psi(alpha) = alpha'*K*alpha + likfun(f), where 
% f = K*alpha+m, and likfun(f) = feval(lik{:},hyp.lik,y,  f,  [],inf).
function [psi,dpsi,f,alpha,dlp,W] = Psi(alpha,m,mvmK,likfun)
  f = mvmK(alpha)+m;
  [lp,dlp,d2lp] = likfun(f); W = -d2lp;
  psi = alpha'*(f-m)/2 - sum(lp);
  if nargout>1, dpsi = mvmK(alpha-dlp); end

% Run IRLS Newton algorithm to optimise Psi(alpha).
function alpha = irls(alpha, m,K,likfun, opt)
  if isfield(opt,'irls_maxit'), maxit = opt.irls_maxit; % max no of Newton steps
  else maxit = 20; end                                           % default value
  if isfield(opt,'irls_Wmin'),  Wmin = opt.irls_Wmin; % min likelihood curvature
  else Wmin = 0.0; end                                           % default value
  if isfield(opt,'irls_tol'),   tol = opt.irls_tol;     % stop Newton iterations
  else tol = 1e-6; end                                           % default value

  smin_line = 0; smax_line = 2;           % min/max line search steps size range
  nmax_line = 10;                          % maximum number of line search steps
  thr_line = 1e-4;                                       % line search threshold
  Psi_line = @(s,alpha,dalpha) Psi(alpha+s*dalpha, m,K.mvm,likfun);% line search
  pars_line = {smin_line,smax_line,nmax_line,thr_line};  % line seach parameters
  search_line = @(alpha,dalpha) brentmin(pars_line{:},Psi_line,5,alpha,dalpha);
  
  f = K.mvm(alpha)+m; [lp,dlp,d2lp] = likfun(f); W = -d2lp; n = size(K,1);
  Psi_new = Psi(alpha,m,K.mvm,likfun);
  Psi_old = Inf;  % make sure while loop starts by the largest old objective val
  it = 0;                          % this happens for the Student's t likelihood
  while Psi_old - Psi_new > tol && it<maxit                       % begin Newton
    Psi_old = Psi_new; it = it+1;
    % limit stepsize
    W = max(W,Wmin); % reduce step size by increasing curvature of problematic W
    [ldB2,solveKiW] = K.fun(W); b = W.*(f-m) + dlp;
    dalpha = b - solveKiW(K.mvm(b)) - alpha;    % Newton direction + line search
    [s_line,Psi_new,n_line,dPsi_new,f,alpha,dlp,W] = search_line(alpha,dalpha);
  end                                                  % end Newton's iterations