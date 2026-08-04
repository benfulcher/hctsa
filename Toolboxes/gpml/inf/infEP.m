function [post nlZ dnlZ] = infEP(hyp, mean, cov, lik, x, y, opt)

% Expectation Propagation approximation to the posterior Gaussian Process.
%
% The (sequential) EP algorithm performs a number of sweeps over the data in
% random order, for better performance when cases are ordered according to the
% targets. During a sweep, every site is locally updated to minimise a local
% divergence measure [1,2] parametrised by a scalar s = [0,1].
% For s=1, one minimises KL(p,q), which amounts to moment matching and
% for s=0, one minimises KL(q,p), which constitutes a 2d joint convex problem
% for log-concave likelihoods.
% The iterations are stopped, when the parameters converge.
%
% Sparse approximations based on inducing inputs (FITC/VFE) [3] are possible
% using the apxSparse functionality.
%
% [1] Tom Minka, Divergence measures and message passing, MSR-TR-2005-173, 2004
% [2] Bui, Yan & Turner, A Unifying Framework for Sparse GP Approximation using
%                        Power EP, 2016, https://arxiv.org/abs/1605.07066
% [3] Naish-Guzman & Holden, The Generalized FITC Approximation, NIPS, 2007
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2017-10-18.
%
% See also infMethods.m, cov/apx.m, gp.m.

persistent last_ttau last_tnu              % keep tilde parameters between calls
tol = 1e-4; max_sweep = 10; min_sweep = 2;     % tolerance to stop EP iterations
if nargin<=6, opt = []; end                        % make opt variable available
if isfield(opt,'s'), s = opt.s; else s = 1.0; end                % default is EP

inf = 'infEP';
n = size(x,1);
if isstruct(cov), K = cov;                   % use provided covariance structure
else K = apx(hyp,cov,x,opt); end               % set up covariance approximation
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv

% init and update drivers for Gaussian posterior approximation
c1 = cov{1}; if isa(c1, 'function_handle'), c1 = func2str(c1); end
sparse = strcmp(c1,'apxSparse') || strcmp(c1,'covFITC');
grid   = strcmp(c1,'apxGrid')   || strcmp(c1,'covGrid');
exact = ~grid && ~sparse;
if s~=0 && s~=1, error('Only opt.s=1 (EP) and opt.s=0 (KL) are possible.'),end
if exact                                         % switch the used functionality
  ep_up = @ep_up_ex; ep_init = @ep_init_ex; ep_marg = @ep_marg_ex;%O(n ^{2,3,0})
elseif sparse
  ep_up = @ep_up_sp; ep_init = @ep_init_sp; ep_marg = @ep_marg_sp;%O(nu^{2,3,2})
elseif grid
  error('EP not implemented for grid covariances.')
end

% A note on naming: variables are given short but descriptive names in 
% accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
% and s2 are mean and variance, nu and tau are natural parameters. A leading t
% means tilde, a subscript _ni means "not i" (for cavity parameters), or _n
% for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

ttau = zeros(n,1); tnu = zeros(n,1);           % init to zero if no better guess
[nlZ,post] = ep_init(K,y,ttau,tnu,lik,hyp,m,inf,s,cov,x);     % init post approx
if all(size(last_ttau)==[n,1])         % try the tilde values from previous call
  [last_nlZ,last_post] = ep_init(K,y,last_ttau,last_tnu,lik,hyp,m,inf,s,cov,x);
  if nlZ > last_nlZ                                         % previous is better
    nlZ = last_nlZ; post = last_post; ttau = last_ttau; tnu = last_tnu;
  end
end

nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
  nlZ_old = nlZ; sweep = sweep+1;
  for i = randperm(n)       % iterate EP updates (in random order) over examples
    [mui,sii] = ep_marg(post,i);                       % obtain marginal moments
    tau_ni = 1/sii-ttau(i);            %  first find the cavity distribution ..
    nu_ni = mui/sii-tnu(i);                         % .. params tau_ni and nu_ni
    if s==1                                    % normal infEP minimising KL(p,q)
      % compute the desired derivatives of the indivdual log partition function
      [lZ, dlZ, d2lZ] = feval(lik{:},hyp.lik,y(i),nu_ni/tau_ni,1/tau_ni,inf);
      ttaui =                     -d2lZ  /(1+d2lZ/tau_ni);    % new tilde params
      ttaui = max(ttaui,0);   % enforce positivity i.e. lower bound ttau by zero
      tnui  = ( dlZ - nu_ni/tau_ni*d2lZ )/(1+d2lZ/tau_ni);
    elseif s==0                                       % infKL minimising KL(q,p)
      [mi,svi] = klmin(lik, hyp.lik, y(i), nu_ni,tau_ni);        % KL projection
      ttaui = 1/svi^2-tau_ni;                                 % new tilde params
      ttaui = max(ttaui,0);   % enforce positivity i.e. lower bound ttau by zero
      tnui  = mi/svi^2-nu_ni;
    end
    post = ep_up(post, ttau(i),tnu(i),i,ttaui-ttau(i),tnui-tnu(i));% update post
    ttau(i) = ttaui; tnu(i) = tnui;                % update the tilde parameters
  end
  % recompute since repeated rank-one updates can destroy numerical precision
  [nlZ,post,alpha,tau_n,nu_n,solveKiW,dhyp] = ...
                                    ep_init(K,y,ttau,tnu,lik,hyp,m,inf,s,cov,x);
  post.L = @(r)-K.P(solveKiW(K.Pt(r))); post.alpha = K.P(alpha);
end

if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
  error('Maximum number of sweeps exceeded in function infEP.')
end
last_ttau = ttau; last_tnu = tnu;                       % remember for next call
if nargout>2                                           % do we want derivatives?
  if s==1                                      % normal infEP minimising KL(p,q)
    dnlZ = dhyp(alpha);                              % covariance-related hypers
    dnlZ.lik = zeros(numel(hyp.lik),1);                        % allocate memory
    for i = 1:numel(hyp.lik)                                 % likelihood hypers
      dlik = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf,i);
      dnlZ.lik(i) = -sum(dlik);
    end
    [junk,dlZ] = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf); % mean hypers
    dnlZ.mean = -dm(dlZ);
  elseif s==0                                         % infKL minimising KL(q,p)
    dnlZ = dhyp(alpha); v = post.diagSigma;          % covariance-related hypers
    dv = -ttau/2;    % at convergence we have df = alpha and dv = -W/2 = -ttau/2
    if sparse
      warning('dnlZ.cov for covSparse and s==0 not yet supported\n')
    else
      [junk,dK] = feval(cov{:}, hyp.cov, x); A = post.A;    % not (yet) scalable
      dnlZ.cov = dnlZ.cov - dK( diag(dv)*A'*(A-A') );
    end
    dnlZ.lik = zeros(numel(hyp.lik),1);                        % allocate memory
    for i = 1:numel(hyp.lik)                                 % likelihood hypers
      dnlZ.lik(i) = -sum( likKL(v,lik,hyp.lik,y,K.mvm(alpha)+m,[],[],i) );
    end
    dnlZ.mean = -dm(alpha);                                        % mean hypers
  end
end
post = struct('alpha',post.alpha, 'L',post.L, 'sW',sqrt(ttau));% clean posterior


function [post,alpha,tau_n,nu_n,ldB2,solveKiW,varargout] = ...
                                                    ep_cavity(K,m,ttau,tnu,post)
  varargout = cell(nargout-6,1);                    % return as much as required
  [ldB2,solveKiW,dW,varargout{:}] = K.fun(ttau);      % functions depending on W
  post.diagSigma = 2*dW;
  alpha = tnu-solveKiW(K.mvm(tnu)+m); post.mu = K.mvm(alpha)+m;
  tau_n = 1./post.diagSigma-ttau; nu_n = post.mu./post.diagSigma-tnu; % cavities

function nlZ = ep_Z(post,alpha,tau_n,nu_n,ldB2,solveKiW,triB,...
                                                       K,m,y,ttau,tnu,lik,hyp,s)
  if s==1                                      % normal infEP minimising KL(p,q)
    lZ = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, 'infEP');
    p = tnu-m.*ttau; q = nu_n-m.*tau_n; r = K.mvm(p);        % auxiliary vectors
    nlZ = ldB2-sum(lZ)-p'*r/2 +r'*solveKiW(r)/2 +(post.diagSigma'*p.^2)/2 ...
       -q'*((ttau./tau_n.*q-2*p).*post.diagSigma)/2-sum(log(1+ttau./tau_n))/2;
  else                                                % infKL minimising KL(q,p)
    lp = likKL(post.diagSigma,lik,hyp.lik,y,post.mu);      % evaluate likelihood
    nlZ = ldB2 -sum(lp) + (alpha'*(post.mu-m)+triB-numel(m))/2;    % upper bound
  end


% Sparse EP using FITC
% refresh the representation of the posterior from initial and site parameters
% to prevent possible loss of numerical precision after many updates
% effort is O(n*nu^2) provided that nu<n
function [nlZ,post,alpha,tau_n,nu_n,solveKiW,dhyp] = ...
                                  ep_init_sp(K,y,ttau,tnu,lik,hyp,m,inf,s,cov,x)
  xu = cov{3}; nu = size(xu,1);
  chol_inv = @(A) rot90(rot90(chol(rot90(rot90(A)))'))\eye(nu);   % chol(inv(A))
  Kuu   = feval(cov{2}{:}, hyp.cov, xu);               % get the building blocks
  Ku = feval(cov{2}{:}, hyp.cov, xu, x);
  diagK = feval(cov{2}{:}, hyp.cov, x, 'diag');
  snu2 = 1e-6*(trace(Kuu)/nu);                 % stabilise by 0.1% of signal std
  R0 = chol_inv(Kuu+snu2*eye(nu));         % initial R, used for refresh O(nu^3)
  V = R0*Ku; d0 = diagK-sum(V.*V,1)';  % initial d, needed for refresh O(n*nu^2)
  post.m = m; post.d0 = d0; post.Ku = Ku; post.R0 = R0;

  if max(norm(tnu),norm(ttau))<1e-10 && nargout<3
    post.diagSigma = diagK;
    nlZ = -sum(feval(lik{:}, hyp.lik, y, m, post.diagSigma, inf));
  else
    [post,alpha,tau_n,nu_n,ldB2,solveKiW,dhyp,L,triB] = ...
                                                   ep_cavity(K,m,ttau,tnu,post);
    nlZ = ep_Z(post,alpha,tau_n,nu_n,ldB2,solveKiW,triB,...
                                                      K,m,y,ttau,tnu,lik,hyp,s);
  end
  if nargout>1
    t  = 1./(1+d0.*ttau);                              % temporary variable O(n)
    d  = d0.*t;                                                           % O(n)
    P  = repmat(t',nu,1).*Ku;                                          % O(n*nu)
    T  = repmat((ttau.*t)',nu,1).*V;              % temporary variable O(n*nu^2)
    R  = chol_inv(eye(nu)+V*T');                                     % O(n*nu^3)
    nn = d.*tnu;                                                          % O(n)
    gg = R0'*(R'*(R*(V*(t.*tnu))));                                    % O(n*nu)
    post.d = d; post.P = P; post.R = R; post.nn = nn; post.gg = gg;
  end

% Sparse EP using FITC
% update the representation of the posterior to reflect modification of the site 
% parameters, effort is O(nu^2)
% old site parameters for site i=1..n: ttaui,        tnui
% new site parameters for site i=1..n: ttaui+dttaui, tnui+dtnui
function post = ep_up_sp(post, ttaui,tnui,i,dttaui,dtnui)
  hi = post.nn(i)+post.m(i) + post.P(:,i)'*post.gg;  % post mean of site i O(nu)
  t = 1+dttaui*post.d(i);
  post.d(i) = post.d(i)/t;                                                % O(1)
  post.nn(i) = post.d(i)*(tnui+dtnui);                                    % O(1)
  r = 1+post.d0(i)*ttaui;
  r = r*r/dttaui + r*post.d0(i);
  v = post.R*(post.R0*post.Ku(:,i));
  r = 1/(r+v'*v);
  if r>0
    post.R = cholupdate(post.R,sqrt( r)*post.R'*v,'-');
  else
    post.R = cholupdate(post.R,sqrt(-r)*post.R'*v,'+');
  end
  T = post.R0'*(post.R'*(post.R*(post.R0*post.P(:,i))));
  post.gg = post.gg + ((dtnui-dttaui*(hi-post.m(i)))/t)*T;             % O(nu^2)
  post.P(:,i) = post.P(:,i)/t;                                           % O(nu)

% Sparse EP using FITC
% O(nu^2) function to compute marginal i of the Gaussian posterior
function [mui,sii] = ep_marg_sp(post,i)
  pi = post.P(:,i); t = post.R*(post.R0*pi);               % temporary variables
  sii = post.d(i) + t'*t; mui = post.nn(i) + pi'*post.gg;    % posterior moments


% Dense EP
% O(n^3) function to compute the parameters of the Gaussian posterior
% approximation, Sigma and mu, and the negative log marginal likelihood, nlZ,
% from the current site parameters, ttau and tnu.
function [nlZ,post,alpha,tau_n,nu_n,solveKiW,dhyp] = ...
                                  ep_init_ex(K,y,ttau,tnu,lik,hyp,m,inf,s,cov,x)
  n = numel(y); post = [];                                    % number of inputs
  if max(norm(tnu),norm(ttau))<1e-10 && nargout<3
    post.mu = m; post.Sigma = K.mvm(eye(n)); post.diagSigma = diag(post.Sigma);
    nlZ = -sum(feval(lik{:}, hyp.lik, y, m, post.diagSigma, inf));
  else
    [post,alpha,tau_n,nu_n,ldB2,solveKiW,dhyp,L,triB] = ...
                                                   ep_cavity(K,m,ttau,tnu,post);
    nlZ = ep_Z(post,alpha,tau_n,nu_n,ldB2,solveKiW,triB,...
                                                      K,m,y,ttau,tnu,lik,hyp,s);
    Kd = K.mvm(eye(n)); post.A = eye(n)-solveKiW(Kd); post.Sigma = Kd*post.A;
  end

% Dense EP
% O(n^2) function to compute the parameters of the Gaussian posterior
% approximation, Sigma and mu, from the current site parameters, 
% old site parameters for site i=1..n: ttaui,        tnui
% new site parameters for site i=1..n: ttaui+dttaui, tnui+dtnui
function post = ep_up_ex(post, ttaui,tnui,i,dttaui,dtnui)
  si = post.Sigma(:,i); sii = si(i); ci = dttaui/(1+dttaui*sii);
  mui = post.mu(i);
  post.Sigma = post.Sigma - ci*(si*si');      % rank-1 update Sigma takes O(n^2)
  post.mu = post.mu - (ci*(mui+sii*dtnui)-dtnui)*si;       % .. and recompute mu

% Dense EP
% O(1) function to obtain the marginal i of the Gaussian posterior
function [mui,sii] = ep_marg_ex(post,i)
  sii = post.Sigma(i,i); mui = post.mu(i);


% KL mode
% Compute the Gaussian Q(f)=N(f|m,s^2) minimising the KL divergence
% KL(Q||P) where P is the product of the cavity distribution q_n(f) and the
% likelihood p(y|f) such that P(f) = 1/Z * q_n(f)*p(y|f).
% The cavity distribution q_n(f) is an unnormalised Gaussian with precision
% parameter tau_n and location parameter nu_n, hence the cavity can be written
% as q_n(f) = exp(nu_n*f-tau_n/2*f^2).
% The minimisation is convex iff. the likelihood p(y|f) is log-concave. The
% optimisation is performed using Newton descent with backtracking line search.
function [m,s,kl] = klmin(lik, hyp, y, nu_n, tau_n)
  ep = 1e-9;                                  % tiny Hessian ridge for stability
  gthresh = 1e-8;                               % gradient convergence threshold
  lik_str = lik{1}; if ~ischar(lik_str), lik_str = func2str(lik_str); end
  if strcmp(lik_str,'likGauss')              % likGauss can be done analytically
    sn2 = exp(2*hyp);
    s = 1/sqrt(1/sn2+tau_n); m = s^2*(nu_n+y/sn2);   % init variables to minimum
  else
    s = 1/sqrt(tau_n); m = nu_n/tau_n;   % init variables to cavity distribution
  end
  ll = likKL(s^2,lik,hyp,y,m);                             % evaluate likelihood
  kl = (s^2+m^2)*tau_n/2 - log(s) - nu_n*m - ll; % init the KL div up to a const
  for i=1:20
    [ll,dm,d2m,dv,d2v,dmdv] = likKL(s^2,lik,hyp,y,m);      % evaluate likelihood
    klold = kl; mold = m; sold = s;                        % remember last value
    kl = (s^2+m^2)*tau_n/2 - log(s) - nu_n*m - ll;   % KL-divergence up to const
    dm = tau_n*m-nu_n-dm;    d2m = tau_n-d2m;           % compute kl derivatives
    ds = s*tau_n-1/s-2*s*dv; d2s = tau_n+1/s^2-2*dv-(2*s)^2*d2v; dmds=-2*s*dmdv;
    detH = ((d2m+ep)*(d2s+ep)-dmds^2);                     % Hessian determinant
    m = m-(dm*(d2s+ep)-ds*dmds)/detH; s = s-(ds*(d2m+ep)-dm*dmds)/detH; % Newton
    for j=1:10                                        % backtracking line search
      if klold>kl, break, end           % we did descend so no need to step back
      m = (m+mold)/2; s = (s+sold)/2;
      kl = (s^2+m^2)*tau_n/2 - log(s) - nu_n*m - likKL(s^2,lik,hyp,y,m);
    end
    d = abs(dm)+abs(dv);                                      % overall gradient
    if j==10, m = mold; s = sold; d = 0; end
    if d<gthresh, break, end
  end

% KL mode
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
    N = 20;                                        % number of quadrature points
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