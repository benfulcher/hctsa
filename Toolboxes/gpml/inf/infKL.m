function [post nlZ dnlZ] = infKL(hyp, mean, cov, lik, x, y)

% Approximation to the posterior Gaussian Process by minimization of the 
% KL-divergence. The function is structurally very similar to infEP; the
% only difference being the local divergence measure minimised.
% In infEP, one minimises KL(p,q) whereas in infKL one minimises KL(q,p)
% locally be iterating over the sites. For log-concave likelihoods, the
% latter minimisation constitutes a 2d joint convex problem, so convergence
% is guaranteed.
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m. See also infMethods.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.
%
% See also INFMETHODS.M.

out = false;
persistent last_ttau last_tnu              % keep tilde parameters between calls
tol = 1e-3; max_sweep = 15; min_sweep = 2;     % tolerance to stop KL iterations

n = size(x,1);
K = feval(cov{:}, hyp.cov, x);                  % evaluate the covariance matrix
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector

% A note on naming: variables are given short but descriptive names in 
% accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
% and s2 are mean and variance, nu and tau are natural parameters. A leading t
% means tilde, a subscript _ni means "not i" (for cavity parameters), or _n
% for a vector of cavity parameters. N(f|mu,Sigma) represents the posterior.

% marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
nlZ0 = -sum(likKL(diag(K),lik,hyp.lik,y,m));
if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
  ttau = zeros(n,1);             % initialize to zero if we have no better guess
  tnu  = zeros(n,1);
  Sigma = K;                     % initialize Sigma and mu, the parameters of ..
  mu = m;                              % .. the Gaussian posterior approximation
  nlZ = nlZ0;
else
  ttau = last_ttau;                    % try the tilde values from previous call
  tnu  = last_tnu;
  [Sigma,mu,L,alpha,nlZ] = klComputeParams(K,y,ttau,tnu,lik,hyp,m);
  if nlZ > nlZ0                                           % if zero is better ..
    ttau = zeros(n,1);                    % .. then initialize with zero instead
    tnu  = zeros(n,1); 
    Sigma = K;                   % initialize Sigma and mu, the parameters of ..
    mu = m;                            % .. the Gaussian posterior approximation
    nlZ = nlZ0;
  end
end

nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
while (nlZ_old-nlZ > tol && sweep < max_sweep) || sweep<min_sweep
  nlZ_old = nlZ; sweep = sweep+1;
  if out, fprintf('Sweep %d, nlZ=%f\n',sweep,nlZ), end
  for i = randperm(n)       % iterate EP updates (in random order) over examples
    tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
    nu_ni = mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni

    ttau_old = ttau(i); tnu_old = tnu(i);  % find the new tilde params, keep old
    [mi,svi] = klmin(lik, hyp.lik, y(i), nu_ni,tau_ni);          % KL projection
    ttau(i) = 1/svi^2-tau_ni;
    ttau(i) = max(ttau(i),0); % enforce positivity i.e. lower bound ttau by zero
    tnu(i) = mi/svi^2-nu_ni;    

    dtt = ttau(i)-ttau_old; dtn = tnu(i)-tnu_old;      % rank-1 update Sigma ..
    si = Sigma(:,i); ci = dtt/(1+dtt*si(i));
    Sigma = Sigma - ci*si*si';
    mu = mu - (ci*(mu(i)+si(i)*dtn)-dtn)*si;               % .. and recompute mu
  end
  % recompute since repeated rank-one updates can destroy numerical precision
  [Sigma,mu,L,alpha,nlZ,A] = klComputeParams(K,y,ttau,tnu,lik,hyp,m);
end

if sweep == max_sweep && nlZ_old-nlZ > tol
  error('maximum number of sweeps exceeded in function infKL')
end

last_ttau = ttau; last_tnu = tnu;                       % remember for next call
post.alpha = alpha; post.sW = sqrt(ttau); post.L = L;  % return posterior params

if nargout>2                                           % do we want derivatives?
  v = diag(Sigma); [lp,df,d2f,dv] = likKL(v,lik,hyp.lik,y,mu);
  dnlZ = hyp;                                   % allocate space for derivatives
  for j=1:length(hyp.cov)                                    % covariance hypers
    dK = feval(cov{:},hyp.cov,x,[],j); AdK = A*dK;
    z = diag(AdK) + sum(A.*AdK,2) - sum(A'.*AdK,1)';
    dnlZ.cov(j) = alpha'*dK*(alpha/2-df) - z'*dv;
  end
  for j=1:length(hyp.lik)                                    % likelihood hypers
    lp_dhyp = likKL(v,lik,hyp.lik,y,K*post.alpha+m,[],[],j);
    dnlZ.lik(j) = -sum(lp_dhyp);
  end
  for j=1:length(hyp.mean)                                         % mean hypers
    dm = feval(mean{:}, hyp.mean, x, j);
    dnlZ.mean(j) = -alpha'*dm;
  end
end

% function to compute the parameters of the Gaussian approximation, Sigma and
% mu, from the current site parameters, ttau and tnu. Also returns L and upper
% bound on negative marginal likelihood.
function [Sigma,mu,L,alpha,nlZ,A] = klComputeParams(K,y,ttau,tnu,lik,hyp,m)
  n = length(tnu);                                    % number of training cases
  sW = sqrt(ttau);                                        % compute Sigma and mu
  L = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
  V = L'\(repmat(sW,1,n).*K);
  Sigma = K - V'*V;
  alpha = tnu-sW.*solve_chol(L,sW.*(K*tnu+m));
  mu = K*alpha+m; v = diag(Sigma);

  A = (eye(n)+K*diag(ttau))\eye(n);                           % A = Sigma*inv(K)
  lp = likKL(v,lik,hyp.lik,y,mu);                          % evaluate likelihood
  nlZ = -sum(lp) - (logdet(A)-alpha'*(mu-m)-trace(A)+n)/2;  % upper bound on -lZ

% We compute the Gaussian Q(f)=N(f|m,s^2) minimising the KL divergence
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
 
% log(det(A)) for det(A)>0 using the LU decomposition of A
function y = logdet(A)
  [L,U] = lu(A); u = diag(U); 
  if prod(sign(u))~=det(L), error('det(A)<=0'), end
  y = sum(log(abs(u)));

% Gaussian smoothed likelihood function; instead of p(y|f)=lik(..,f,..) compute
%   log likKL(f) = int log lik(..,t,..) N(f|t,v), where
%     v   .. marginal variance = (positive) smoothing width, and
%     lik .. lik function such that feval(lik{:},varargin{:}) yields a result.
% All return values are separately integrated using Gaussian-Hermite quadrature.
% Gaussian smoothed likelihood function; instead of p(y|f)=lik(..,f,..) compute
%   likKL(f) = int lik(..,t,..) N(f|t,v), where
%     v   .. marginal variance = (positive) smoothing width, and
%     lik .. lik function such that feval(lik{:},varargin{:}) yields a result.
% All return values are separately integrated using Gaussian-Hermite quadrature.
function [ll,df,d2f,dv,d2v,dfdv] = likKL(v, lik, varargin)
  N = 20;                                          % number of quadrature points
  [t,w] = gauher(N);      % location and weights for Gaussian-Hermite quadrature
  f = varargin{3};                               % obtain location of evaluation
  sv = sqrt(v);                                                % smoothing width
  ll = 0; df = 0; d2f = 0; dv = 0; d2v = 0; dfdv = 0;    % init return arguments
  for i=1:N                                            % use Gaussian quadrature
    varargin{3} = f + sv*t(i);   % coordinate transform of the quadrature points
    [lp,dlp,d2lp] = feval(lik{:},varargin{1:3},[],'infLaplace',varargin{6:end});
    if nargout>0,     ll  = ll  + w(i)*lp;               % value of the integral
      if nargout>1,   df  = df  + w(i)*dlp;             % derivative w.r.t. mean
        if nargout>2, d2f = d2f + w(i)*d2lp;        % 2nd derivative w.r.t. mean
          if nargout>3                              % derivative w.r.t. variance
            ai = t(i)./(2*sv+eps); dvi = dlp.*ai; dv = dv + w(i)*dvi; % no 0 div
            if nargout>4                        % 2nd derivative w.r.t. variance
              d2v = d2v + w(i)*(d2lp.*(t(i)^2/2)-dvi)./(v+eps)/2;     % no 0 div
              if nargout>5                            % mixed second derivatives
                dfdv = dfdv + w(i)*(ai.*d2lp);
              end
            end
          end
        end
      end
    end
  end
  % likLaplace can be done analytically
  lik_str = lik{1}; if ~ischar(lik_str), lik_str = func2str(lik_str); end
  if strcmp(lik_str,'likLaplace')
    b = exp(varargin{1})/sqrt(2); y = varargin{2}; sv = sqrt(v);
    mu = (f-y)/b; z = (f-y)./sv;
    Nz = exp(-z.^2/2)/sqrt(2*pi);
    Cz = (1+erf(z/sqrt(2)))/2;
    ll = (1-2*Cz).*mu - 2/b*sv.*Nz - log(2*b);
    df = (1-2*Cz)/b;
    d2f = -2*Nz./(b*(sv+eps));
    dv = d2f/2;
    d2v = (z.*z-1)./(v+eps).*d2f/4;
    dfdv = -z.*d2f./(2*sv+eps);
  end
