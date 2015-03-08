function [post nlZ dnlZ] = infEP(hyp, mean, cov, lik, x, y)

% Expectation Propagation approximation to the posterior Gaussian Process.
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m. See also infMethods.m. In the EP algorithm, the sites are 
% updated in random order, for better performance when cases are ordered
% according to the targets.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-09-13.
%
% See also INFMETHODS.M.

persistent last_ttau last_tnu              % keep tilde parameters between calls
tol = 1e-4; max_sweep = 10; min_sweep = 2;     % tolerance to stop EP iterations

inf = 'infEP';
n = size(x,1);
K = feval(cov{:}, hyp.cov, x);                  % evaluate the covariance matrix
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector

% A note on naming: variables are given short but descriptive names in 
% accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
% and s2 are mean and variance, nu and tau are natural parameters. A leading t
% means tilde, a subscript _ni means "not i" (for cavity parameters), or _n
% for a vector of cavity parameters. N(f|mu,Sigma) is the posterior.

% marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
nlZ0 = -sum(feval(lik{:}, hyp.lik, y, m, diag(K), inf));
if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
  ttau = zeros(n,1); tnu  = zeros(n,1);        % init to zero if no better guess
  Sigma = K;                     % initialize Sigma and mu, the parameters of ..
  mu = m; nlZ = nlZ0;                  % .. the Gaussian posterior approximation
else
  ttau = last_ttau; tnu  = last_tnu;   % try the tilde values from previous call
  [Sigma,mu,L,alpha,nlZ] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
  if nlZ > nlZ0                                           % if zero is better ..
    ttau = zeros(n,1); tnu  = zeros(n,1);       % .. then init with zero instead
    Sigma = K;                   % initialize Sigma and mu, the parameters of ..
    mu = m; nlZ = nlZ0;                % .. the Gaussian posterior approximation
  end
end

nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
  nlZ_old = nlZ; sweep = sweep+1;
  for i = randperm(n)       % iterate EP updates (in random order) over examples
    tau_ni = 1/Sigma(i,i)-ttau(i);      %  first find the cavity distribution ..
    nu_ni = mu(i)/Sigma(i,i)-tnu(i);                % .. params tau_ni and nu_ni

    % compute the desired derivatives of the indivdual log partition function
    [lZ, dlZ, d2lZ] = feval(lik{:}, hyp.lik, y(i), nu_ni/tau_ni, 1/tau_ni, inf);
    ttau_old = ttau(i); tnu_old = tnu(i);  % find the new tilde params, keep old
    ttau(i) =                     -d2lZ  /(1+d2lZ/tau_ni);
    ttau(i) = max(ttau(i),0); % enforce positivity i.e. lower bound ttau by zero
    tnu(i)  = ( dlZ - nu_ni/tau_ni*d2lZ )/(1+d2lZ/tau_ni);

    dtt = ttau(i)-ttau_old; dtn = tnu(i)-tnu_old;      % rank-1 update Sigma ..
    si = Sigma(:,i); ci = dtt/(1+dtt*si(i));
    Sigma = Sigma - ci*si*si';                         % takes 70% of total time
    mu = mu - (ci*(mu(i)+si(i)*dtn)-dtn)*si;               % .. and recompute mu
  end
  % recompute since repeated rank-one updates can destroy numerical precision
  [Sigma,mu,L,alpha,nlZ] = epComputeParams(K, y, ttau, tnu, lik, hyp, m, inf);
end

if sweep == max_sweep && abs(nlZ-nlZ_old) > tol
  error('maximum number of sweeps exceeded in function infEP')
end

last_ttau = ttau; last_tnu = tnu;                       % remember for next call
post.alpha = alpha; post.sW = sqrt(ttau); post.L = L;  % return posterior params

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
  nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
  sW = sqrt(ttau);
  F = alpha*alpha'-repmat(sW,1,n).*solve_chol(L,diag(sW));   % covariance hypers
  for i=1:length(hyp.cov)
    dK = feval(cov{:}, hyp.cov, x, [], i);
    dnlZ.cov(i) = -sum(sum(F.*dK))/2;
  end
  for i = 1:numel(hyp.lik)                                   % likelihood hypers
    dlik = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf, i);
    dnlZ.lik(i) = -sum(dlik);
  end
  [junk,dlZ] = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);% mean hyps
  for i = 1:numel(hyp.mean)
    dm = feval(mean{:}, hyp.mean, x, i);
    dnlZ.mean(i) = -dlZ'*dm;
  end
end

% function to compute the parameters of the Gaussian approximation, Sigma and
% mu, and the negative log marginal likelihood, nlZ, from the current site
% parameters, ttau and tnu. Also returns L (useful for predictions).
function [Sigma,mu,L,alpha,nlZ] = epComputeParams(K,y,ttau,tnu,lik,hyp,m,inf)
  n = length(y);                                      % number of training cases
  sW = sqrt(ttau);                                        % compute Sigma and mu
  L = chol(eye(n)+sW*sW'.*K);                            % L'*L=B=eye(n)+sW*K*sW
  V = L'\(repmat(sW,1,n).*K);
  Sigma = K - V'*V;
  alpha = tnu-sW.*solve_chol(L,sW.*(K*tnu+m));
  mu = K*alpha+m; v = diag(Sigma);

  tau_n = 1./diag(Sigma)-ttau;             % compute the log marginal likelihood
  nu_n  = mu./diag(Sigma)-tnu;                    % vectors of cavity parameters
  lZ = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);
  p = tnu-m.*ttau; q = nu_n-m.*tau_n;                        % auxiliary vectors
  nlZ = sum(log(diag(L))) - sum(lZ) - p'*Sigma*p/2 + (v'*p.^2)/2 ...
      - q'*((ttau./tau_n.*q-2*p).*v)/2 - sum(log(1+ttau./tau_n))/2;
