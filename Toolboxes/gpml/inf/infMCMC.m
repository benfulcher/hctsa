function [post nlZ dnlZ] = infMCMC(hyp, mean, cov, lik, x, y, par)

% Markov Chain Monte Carlo (MCMC) sampling from posterior and 
% Annealed Importance Sampling (AIS) for marginal likelihood estimation.
%
% The algorithms are not to be used as a black box, since the acceptance rate
% of the samplers need to be carefully monitored. Also, there are no derivatives 
% of the marginal likelihood available.
%
% There are additional parameters:
%   - par.sampler   switch between the samplers 'hmc' or 'ess'
%   - par.Nsample   num of samples
%   - par.Nskip     num of steps out of which one sample kept
%   - par.Nburnin   num of burn in samples (corresponds to Nskip*Nburning steps)
%   - par.Nais      num of AIS runs to remove finite temperature bias
%   - par.st        spherical Gaussian width of the posterior KDE
% Default values are sampler=hmc, Nsample=200, Nskip=40, Nburnin=10, Nais=3,
%                                 st=0.001*sqrt(sum(diag(K))/n).
%
% The Hybrid Monte Carlo Sampler (HMC) is implemented as described in the
% technical report: Probabilistic Inference using MCMC Methods by Radford Neal,
% CRG-TR-93-1, 1993.
%
%     Instead of sampling from f ~ 1/Zf * N(f|m,K) P(y|f), we use a 
%     parametrisation in terms of alpha = inv(K)*(f-m) and sample from 
%     alpha ~ P(a) = 1/Za * N(a|0,inv(K)) P(y|K*a+m) to increase numerical 
%     stability since log P(a) = -(a'*K*a)/2 + log P(y|K*a+m) + C and its 
%     gradient can be computed safely.
%
%     The leapfrog stepsize as a time discretisation stepsize comes with a 
%     tradeoff:
%       - too small: frequently accept, slow exploration, accurate dynamics
%       - too large: seldomly reject, fast exploration, inaccurate dynamics
%     In order to balance between the two extremes, we adaptively adjust the
%     stepsize in order to keep the acceptance rate close to a target acceptance
%     rate. Taken from http://deeplearning.net/tutorial/hmc.html.
%
%     The code issues a warning in case the overall acceptance rate did deviate
%     strongly from the target. This can indicate problems with the sampler.
%
% The Elliptical Slice Sampler (ESS) is a straight implementation that is
% inspired by the code shipped with the paper: Elliptical slice sampling by
% Iain Murray, Ryan Prescott Adams and David J.C. MacKay, AISTATS 2010.
%
% Annealed Importance Sampling (AIS) to determine the marginal likelihood is
% described in the technical report: AIS, Radford Neal, 1998.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-05-09.
%
% See also infMethods.m.

if nargin<7, par = []; end                         % analyse parameter structure
if isfield(par,'sampler'), alg=par.sampler; else alg='hmc';   end
if isfield(par,'Nsample'), N  =par.Nsample; else N  =200;     end
if isfield(par,'Nskip'),   Ns =par.Nskip;   else Ns = 40;     end
if isfield(par,'Nburnin'), Nb =par.Nburnin; else Nb = 10;     end
if isfield(par,'Nais'),    R  =par.Nais;    else R  =  3;     end

K = feval(cov{:}, hyp.cov, x);                  % evaluate the covariance matrix
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector
n = size(K,1); st_def = 0.001*sqrt(sum(diag(K))/n);              % default value
if isfield(par,'st'),      st =par.st;      else st = st_def; end
[cK,fail] = chol(K+st^2*eye(n));        % try an ordinary Cholesky decomposition
if fail, error('Try increasing the par.st variable.'), end
T = (N+Nb)*Ns;                                         % overall number of steps
[alpha,Na] = sample(K,cK,m,y,lik,hyp.lik, N,Nb,Ns, alg);  % sample w/o annealing
post.alpha = alpha; al = sum(alpha,2)/N;
post.L  = -(cK\(cK'\eye(n)) + al*al' - alpha*alpha'/N);    % inv(K) - cov(alpha)
post.sW = [];
post.acceptance_rate_MCMC = Na/T;                  % additional output parameter
if nargout>1                                      % annealed importance sampling
  % discrete time t from 1 to T and temperature tau from tau(1)=0 to tau(T)=1
  taus = [zeros(1,Nb),linspace(0,1,N)].^4; % annealing schedule, effort at start
  % the fourth power idea is taken from the Kuss&Rasmussen paper, 2005 JMLR
  % Z(t) := \int N(f|m,K) lik(f)^taus(t) df hence Z(1) = 1 and Z(T) = Z
  % Z = Z(T)/Z(1) = prod_t Z(t)/Z(t-1); 
  % ln Z(t)/Z(t-1) = ( tau(t)-tau(t-1) ) * loglik(f_t)
  lZ = zeros(R,1); dtaus = diff(taus(Nb+(1:N)));        % we have: sum(dtaus)==1
  for r=1:R
    [A,Na] = sample(K,cK,m,y,lik,hyp.lik, N,Nb,Ns, alg, taus);             % AIS
    for t=2:N                              % evaluate the likelihood sample-wise
      lp = feval(lik{:},hyp.lik,y,K*A(:,t)+m,[],'infLaplace');
      lZ(r) = lZ(r)+dtaus(t-1)*sum(lp);
    end
    post.acceptance_rate_AIS(r) = Na/T;           % additional output parameters
  end
  nlZ = log(R)-logsumexp(lZ);  % remove finite temperature bias, softmax average
  if nargout>2                % marginal likelihood derivatives are not computed
    dnlZ = struct('cov',0*hyp.cov, 'mean',0*hyp.mean, 'lik',0*hyp.lik);
  end
end

%% choose between HMC and ESS depending on the alg string
function [alpha,Na] = sample(K,cK,m,y,lik,hyp, N,Nb,Ns, alg, varargin)
  if strcmp(lik,'likGauss')
    [alpha,Na] = sample_gauss(K, m,y,hyp, N,Nb,Ns, varargin{:});
  else
    if strcmpi(alg,'hmc')
      [alpha,Na] = sample_hmc(K, m,y,lik,hyp, N,Nb,Ns, varargin{:});
    else
      [alpha,Na] = sample_ess(cK,m,y,lik,hyp, N,Nb,Ns, varargin{:});
    end
  end

%% sample from Gaussian posterior
function [alpha,Na] = sample_gauss(K, m,y,hyp, N,Nb,Ns, taus)
  Na = (N+Nb)*Ns; sn2 = exp(2*hyp); n = size(K,1);
  if nargin<8
    S = K + K*K/sn2; alpha = chol(S)\randn(n,N) + repmat(S\(K*(y-m)/sn2),1,N);
  else
    alpha = zeros(n,N);
    for t=1:N
      St = K + K*K*(taus(Nb+t)/sn2);
      alpha(:,t) = chol(St)\randn(n,1) + St\(K*(y-m)*taus(Nb+t)/sn2);
    end
  end

%% sample using elliptical slices
function [alpha,T] = sample_ess(cK,m,y,lik,hyp, N,Nb,Ns, taus)
  if nargin>=9, tau = taus(1); else tau = 1; end       % default is no annealing
  T = (N+Nb)*Ns;                                       % overall number of steps
  F = zeros(size(m,1),N);
  for t=1:T
    if nargin>=9, tau = taus(1+floor((t-1)/Ns)); end   % parameter from schedule
    if t==1, f=0*m; l=sum(feval(lik{:},hyp,y,f+m)); end  % init sample f & lik l
    r = cK'*randn(size(m));                          % random sample from N(0,K)
    [f,l] = sample_ess_step(f,l,r,m,y,lik,hyp,tau);
    if mod(t,Ns)==0                                   % keep one state out of Ns
      if t/Ns>Nb, F(:,t/Ns-Nb) = f; end              % wait for Nb burn-in steps
    end
  end
  alpha = cK\(cK'\F);

%% elliptical slice sampling: one step
function [f,l] = sample_ess_step(f,l,r,m,y,lik,hyp,tau)
  if nargin<8, tau=1; end, if tau>1, tau=1; end, if tau<0, tau=0; end
  h = log(rand) + tau*l;                                  % acceptance threshold
  a = rand*2*pi; amin = a-2*pi; amax = a;                % bracket whole ellipse
  k = 0;                                                       % emergency break
  while true % slice sampling loop; f for proposed angle diff; check if on slice
    fp = f*cos(a) + r*sin(a);                    % move on ellipsis defined by r
    l = sum(feval(lik{:},hyp,y,fp+m));
    if tau*l>h || k>20, break, end  % exit if new point is on slice or exhausted
    if a>0, amax=a; elseif a<0, amin=a; end     % shrink slice to rejected point
    a = rand*(amax-amin) + amin; k = k+1;  % propose new angle difference; break
  end
  f = fp;                                                               % accept

%% sample using Hamiltonian dynamics as proposal algorithm
function [alpha,Na] = sample_hmc(K,m,y,lik,hyp, N,Nb,Ns, taus)
  % use adaptive stepsize rule to enforce a specific acceptance rate 
  epmin = 1e-6; % minimum leapfrog stepsize
  epmax = 9e-1; % maximum leapfrog stepsize
  ep    = 1e-2; % initial leapfrog stepsize
  acc_t = 0.9;  % target acceptance rate
  acc   = 0;    % current acceptance rate
  epinc = 1.02; % increase factor of stepsize if acceptance rate is below target
  epdec = 0.98; % decrease factor of stepsize if acceptance rate is above target
  lam   = 0.01; %  exponential moving average computation of the acceptance rate
          % 2/(3*lam) steps half height; lam/nstep = 0.02/33, 0.01/66, 0.005/133
  l = 20;       % number of leapfrog steps to perform for one step
  n = size(K,1);
  T = (N+Nb)*Ns;                                       % overall number of steps
  alpha = zeros(n,N);                                            % sample points
  al = zeros(n,1);                                            % current position
  if nargin>=9, tau = taus(1); else tau = 1; end       % default is no annealing
  [gold,eold] = E(al,K,m,y,lik,hyp,tau);              % initial energy, gradient
  Na = 0;                                            % number of accepted points
  for t=1:T
    if nargin>=9, tau = taus(1+floor((t-1)/Ns)); end   % parameter from schedule
    p = randn(n,1);                                    % random initial momentum
    q = al;
    g = gold;
    Hold = (p'*p)/2 + eold;                     % H = Ekin + Epot => Hamiltonian
    for ll=1:l                       % leapfrog discretization steps, Euler like
      p = p - (ep/2)*g;              % half step in momentum p
      q = q + ep*p;                  % full step in position q
      g = E(q,K,m,y,lik,hyp,tau);    % compute new gradient  g
      p = p - (ep/2)*g;              % half step in momentum p
    end
    [g,e] = E(q,K,m,y,lik,hyp,tau);
    H = (p'*p)/2 + e;               % recompute Hamiltonian
    acc = (1-lam)*acc;              % decay current acceptance rate
    if log(rand) < Hold-H           % accept with p = min(1,exp(Hold-H))
      al         = q;               % keep new state,
      gold       = g;               % gradient
      eold       = e;               % and potential energy
      acc = acc + lam;              % increase rate due to acceptance
      Na = Na+1;                    % increase number accepted steps
    end
    if acc>acc_t       % too large acceptance rate => increase leapfrog stepsize
      ep = epinc*ep;
    else               % too small acceptance rate => decrease leapfrog stepsize
      ep = epdec*ep;
    end
    if ep<epmin, ep = epmin; end             % clip stepsize ep to [epmin,epmax]
    if ep>epmax, ep = epmax; end
    if mod(t,Ns)==0                                   % keep one state out of Ns
      if t/Ns>Nb, alpha(:,t/Ns-Nb) = al; end         % wait for Nb burn-in steps
    end
  end
  if Na/T<acc_t*0.9 || 1.07*acc_t<Na/T  % Acceptance rate in the right ballpark?
    fprintf('The acceptance rate %1.2f%% is not within',100*Na/T)
    fprintf(' [%1.1f, %1.1f]%%\n', 100*acc_t*0.9, 100*acc_t*1.07)
    if nargin<9
      warning('Bad (HMC) acceptance rate')
    else
      warning('Bad (AIS) acceptance rate')
    end
  end

%% potential energy function value and its gradient
function [g,e] = E(al,K,m,y,lik,hyp,tau)
  % E(f) = E(al) = al'*K*al/2 - tau*loglik(f), f = K*al+m
  % g(f) = g(al) = al - tau*dloglik(f) => dE/dal = K*( al-dloglik(f) )
  if nargin<7, tau=1; end, if tau>1, tau=1; end, if tau<0, tau=0; end
  Kal = K*al;
  [lp,dlp] = feval(lik{:},hyp,y,Kal+m,[],'infLaplace');
  g = Kal-K*(tau*dlp);
  if nargout>1, e = (al'*Kal)/2 - tau*sum(lp); end

%% y = logsumexp(x) = log(sum(exp(x(:))) avoiding overflow
function y=logsumexp(x)
  mx = max(x(:));
  y = log(sum(exp(x(:)-mx)))+mx;
