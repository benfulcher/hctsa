function [post nlZ dnlZ] = infKL(hyp, mean, cov, lik, x, y, opt)

% Approximation to the posterior Gaussian Process by direct minimization of the 
% KL-divergence.
%
% The optimisation procedure is an implementation of the CVI algorithm for GP
% models, based on Algorithm 1 and Equation 13 in [1].
% To estimate gradients of E[log p(y|f)], we use Gauss-Hermite quadrature.
%
% [1] Emtiyaz Khan and Wu Lin, Conjugate-Computation Variational Inference:
% Converting Variational Inference in Non-Conjugate Models to Inferences in
% Conjugate Models, AISTATS, 2017.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Emtiyaz Khan (RIKEN) and Wu Lin (RIKEN) 2017-08-18
%                     with some changes by Hannes Nickisch 2017-10-22.
%
% See also infMethods.m, cov/apx.m, gp.m.

c1 = cov{1}; if isa(c1, 'function_handle'), c1 = func2str(c1); end
spars = strcmp(c1,'apxSparse') || strcmp(c1,'covFITC');
grid  = strcmp(c1,'apxGrid')   || strcmp(c1,'covGrid');
state = strcmp(c1,'apxState');
exact = ~grid && ~spars && ~state;

if nargin<7, opt = []; end
if isfield(opt,'step_size'), beta = opt.step_size; else beta = 0.1; end
if isfield(opt,'tol'), tol = opt.tol; else tol = 1e-4; end
if isfield(opt,'max_iters'), kmax = opt.max_iters; else kmax = 100; end

n = size(x,1);
if isstruct(cov), K = cov;                   % use provided covariance structure
else K = apx(hyp,cov,x,opt); end               % set up covariance approximation
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv

tlam1 = 1e-3*ones(n,1); tlam2 = 1e-3*ones(n,1);  % init natural params of GP apx

for k = 1:kmax    % iterate using eq 13 and alg 1 in Khan and Lin, AISTATS, 2017
  yp = tlam1./tlam2; sW = sqrt(abs(tlam2));        % Gaussian pseudo-observation
  [ldB2,solveKiW,dW,dhyp,L,triB] = K.fun(sW.*sW);         % GP regression update
  v = 2*dW; alpha = solveKiW(yp-m); mu = K.mvm(alpha)+m; % marginal var and mean
  [ll, df, junk, dv] = likKL(v, lik, hyp.lik, y, mu);               % likelihood
  nlZ = ldB2 -sum(ll) + (alpha'*(mu-m)+triB-numel(m))/2;   % marginal likelihood
  tlam1 = (1-beta).*tlam1 + beta.*( df-2*dv.*mu );    % natural parameter update
  tlam2 = (1-beta).*tlam2 + beta.*(   -2*dv     );
  if k==1; nlZ_old = nlZ; end                             % diagnose convergence
  if abs(nlZ-nlZ_old)<tol && k>1, break, end
  nlZ_old = nlZ;
  fprintf('%03d) nlZ=%1.4f\n',k,nlZ)
end

if k>=kmax, warning('Max number of iterations reached.\n'), end
post.sW = sW; post.L = @(r)-K.P(solveKiW(K.Pt(r))); post.alpha = K.P(alpha);

if nargout>2                                           % do we want derivatives?
  W = sW.*sW; [ldB2,solveKiW,dW,dhyp] = K.fun(W);
  dnlZ = dhyp(alpha); v = 2*dW;                      % covariance-related hypers
  dv = -W/2;                   % at convergence we have df = alpha and dv = -W/2
  if ~exact
    warning('Derivatives not yet scalable for apx*.')
  else
    Kd = K.mvm(eye(n)); A = eye(n)-solveKiW(Kd);            % not (yet) scalable
    [junk,dK] = feval(cov{:}, hyp.cov, x);
    dnlZ.cov = dnlZ.cov - dK( diag(dv)*A'*(A-A') );
  end
  dnlZ.lik = zeros(1,numel(hyp.lik));                          % allocate memory
  for i = 1:numel(hyp.lik)                                   % likelihood hypers
    dnlZ.lik(i) = -sum( likKL(v,lik,hyp.lik,y,K.mvm(alpha)+m,[],[],i) );
  end
  dnlZ.mean = -dm(alpha);                                          % mean hypers
end

% Gaussian smoothed likelihood function; instead of p(y|f)=lik(..,f,..) compute
%   log likKL(f) = int log lik(..,t,..) N(f|t,v) dt, where
%     v   .. marginal variance = (positive) smoothing width, and
%     lik .. lik function such that feval(lik{:},varargin{:}) yields a result.
% All return values are separately integrated using Gaussian-Hermite quadrature.
function [ll,df,d2f,dv,d2v,dfdv] = likKL(v, lik, varargin)
  f = varargin{3};                               % obtain location of evaluation
  sv = sqrt(v);                                                % smoothing width
  lik_str = lik{1}; if ~ischar(lik_str), lik_str = func2str(lik_str); end
  if strcmp(lik_str,'likLaplace')          % likLaplace can be done analytically
    b = exp(varargin{1})/sqrt(2); y = varargin{2};
    mu = (f-y)/b; z = (f-y)./sv;
    Nz = exp(-z.^2/2)/sqrt(2*pi);
    Cz = (1+erf(z/sqrt(2)))/2;
    ll = (1-2*Cz).*mu - 2/b*sv.*Nz - log(2*b);
    df = (1-2*Cz)/b;
    d2f = -2*Nz./(b*(sv+eps));
    dv = d2f/2;
    d2v = (z.*z-1)./(v+eps).*d2f/4;
    dfdv = -z.*d2f./(2*sv+eps);
  else
    N = 50;                                        % number of quadrature points
    [t,w] = gauher(N);    % location and weights for Gaussian-Hermite quadrature
    ll = 0; df = 0; d2f = 0; dv = 0; d2v = 0; dfdv = 0;  % init return arguments
    for i=1:N                                          % use Gaussian quadrature
      varargin{3} = f + sv*t(i); % coordinate transform of the quadrature points
      [lp,dlp,d2lp]=feval(lik{:},varargin{1:3},[],'infLaplace',varargin{6:end});
      ll   = ll  + w(i)*lp;                              % value of the integral
      df   = df  + w(i)*dlp;                              % derivative wrt. mean
      d2f  = d2f + w(i)*d2lp;                         % 2nd derivative wrt. mean
      ai = t(i)./(2*sv+eps); dvi = dlp.*ai; dv = dv+w(i)*dvi;   % deriv wrt. var
      d2v  = d2v + w(i)*(d2lp.*(t(i)^2/2)-dvi)./(v+eps)/2;  % 2nd deriv wrt. var
      dfdv = dfdv + w(i)*(ai.*d2lp);                  % mixed second derivatives
    end
  end