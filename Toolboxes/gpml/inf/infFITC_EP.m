function [post nlZ dnlZ] = infFITC_EP(hyp, mean, cov, lik, x, y)

% FITC-EP approximation to the posterior Gaussian process. The function is
% equivalent to infEP with the covariance function:
%   Kt = Q + G; G = diag(g); g = diag(K-Q);  Q = Ku'*inv(Kuu + snu2*eye(nu))*Ku;
% where Ku and Kuu are covariances w.r.t. to inducing inputs xu and
% snu2 = sn2/1e6 is the noise of the inducing inputs. We fixed the standard
% deviation of the inducing inputs snu to be a one per mil of the measurement 
% noise's standard deviation sn. In case of a likelihood without noise
% parameter sn2, we simply use snu2 = 1e-6.
% For details, see The Generalized FITC Approximation, Andrew Naish-Guzman and
%                  Sean Holden, NIPS, 2007.
%
% The implementation exploits the Woodbury matrix identity
%   inv(Kt) = inv(G) - inv(G)*Ku'*inv(Kuu+Ku*inv(G)*Ku')*Ku*inv(G)
% in order to be applicable to large datasets. The computational complexity
% is O(n nu^2) where n is the number of data points x and nu the number of
% inducing inputs in xu.
% The posterior N(f|h,Sigma) is given by h = m+mu with mu = nn + P'*gg and
% Sigma = inv(inv(K)+diag(W)) = diag(d) + P'*R0'*R'*R*R0*P. Here, we use the
% site parameters: b,w=$b,\pi$=tnu,ttau, P=$P'$, nn=$\nu$, gg=$\gamma$
%             
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with covFITC.
%
% The inducing points can be specified through 1) the 2nd covFITC parameter or
% by 2) providing a hyp.xu hyperparameters. Note that 2) has priority over 1).
% In case 2) is provided and derivatives dnlZ are requested, there will also be
% a dnlZ.xu field allowing to optimise w.r.t. to the inducing points xu. However
% the derivatives dnlZ.xu can only be computed for one of the following eight
% covariance functions: cov{Matern|PP|RQ|SE}{iso|ard}.
%
% Copyright (c) by Hannes Nickisch, 2013-10-29.
%
% See also INFMETHODS.M, COVFITC.M.

persistent last_ttau last_tnu              % keep tilde parameters between calls
tol = 1e-4; max_sweep = 20; min_sweep = 2;     % tolerance to stop EP iterations
inf = 'infEP';
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if ~strcmp(cov1,'covFITC'); error('Only covFITC supported.'), end    % check cov
if isfield(hyp,'xu'), cov{3} = hyp.xu; end  % hyp.xu is provided, replace cov{3}

[diagK,Kuu,Ku] = feval(cov{:}, hyp.cov, x);         % evaluate covariance matrix
if ~isempty(hyp.lik)                          % hard coded inducing inputs noise
  sn2 = exp(2*hyp.lik(end)); snu2 = 1e-6*sn2;               % similar to infFITC
else
  snu2 = 1e-6;
end
[n, D] = size(x); nu = size(Kuu,1);
m = feval(mean{:}, hyp.mean, x);                      % evaluate the mean vector

rot180   = @(A)   rot90(rot90(A));                     % little helper functions
chol_inv = @(A) rot180(chol(rot180(A))')\eye(nu);                 % chol(inv(A))
R0 = chol_inv(Kuu+snu2*eye(nu));           % initial R, used for refresh O(nu^3)
V = R0*Ku; d0 = diagK-sum(V.*V,1)';    % initial d, needed for refresh O(n*nu^2)

% A note on naming: variables are given short but descriptive names in 
% accordance with Rasmussen & Williams "GPs for Machine Learning" (2006): mu
% and s2 are mean and variance, nu and tau are natural parameters. A leading t
% means tilde, a subscript _ni means "not i" (for cavity parameters), or _n
% for a vector of cavity parameters.

% marginal likelihood for ttau = tnu = zeros(n,1); equals n*log(2) for likCum*
nlZ0 = -sum(feval(lik{:}, hyp.lik, y, m, diagK, inf));
if any(size(last_ttau) ~= [n 1])      % find starting point for tilde parameters
  ttau = zeros(n,1);             % initialize to zero if we have no better guess
  tnu  = zeros(n,1);
  [d,P,R,nn,gg] = epfitcRefresh(d0,Ku,R0,V, ttau,tnu); % compute initial repres.
  nlZ = nlZ0;
else
  ttau = last_ttau;                    % try the tilde values from previous call
  tnu  = last_tnu;
  [d,P,R,nn,gg] = epfitcRefresh(d0,Ku,R0,V, ttau,tnu); % compute initial repres.
  nlZ = epfitcZ(d,P,R,nn,gg,ttau,tnu,d0,R0,Ku,y,lik,hyp,m,inf);
  if nlZ > nlZ0                                           % if zero is better ..
    ttau = zeros(n,1);                    % .. then initialize with zero instead
    tnu  = zeros(n,1);
    [d,P,R,nn,gg] = epfitcRefresh(d0,Ku,R0,V, ttau,tnu);       % initial repres.
    nlZ = nlZ0;
  end
end

nlZ_old = Inf; sweep = 0;               % converged, max. sweeps or min. sweeps?
while (abs(nlZ-nlZ_old) > tol && sweep < max_sweep) || sweep<min_sweep
  nlZ_old = nlZ; sweep = sweep+1;
  for i = randperm(n)       % iterate EP updates (in random order) over examples
    pi = P(:,i); t = R*(R0*pi);                            % temporary variables
    sigmai = d(i) + t'*t; mui = nn(i) + pi'*gg;           % post moments O(nu^2)
    
    tau_ni = 1/sigmai-ttau(i);          %  first find the cavity distribution ..
    nu_ni = mui/sigmai+m(i)*tau_ni-tnu(i);          % .. params tau_ni and nu_ni
   
    % compute the desired derivatives of the indivdual log partition function
    [lZ, dlZ, d2lZ] = feval(lik{:}, hyp.lik, y(i), nu_ni/tau_ni, 1/tau_ni, inf);
    ttaui =                            -d2lZ  /(1+d2lZ/tau_ni);
    ttaui = max(ttaui,0);     % enforce positivity i.e. lower bound ttau by zero
    tnui  = ( dlZ + (m(i)-nu_ni/tau_ni)*d2lZ )/(1+d2lZ/tau_ni);
    [d,P(:,i),R,nn,gg,ttau,tnu] = ...                    % update representation
              epfitcUpdate(d,P(:,i),R,nn,gg, ttau,tnu,i,ttaui,tnui, m,d0,Ku,R0);
  end
  % recompute since repeated rank-one updates can destroy numerical precision
  [d,P,R,nn,gg] = epfitcRefresh(d0,Ku,R0,V, ttau,tnu);
  [nlZ,nu_n,tau_n] = epfitcZ(d,P,R,nn,gg,ttau,tnu,d0,R0,Ku,y,lik,hyp,m,inf);
end

if sweep == max_sweep
  warning('maximum number of sweeps reached in function infEP')
end
last_ttau = ttau; last_tnu = tnu;                       % remember for next call

post.sW = sqrt(ttau);                  % unused for FITC_EP prediction with gp.m
dd = 1./(d0+1./ttau);
alpha = tnu./(1+d0.*ttau);
RV = R*V; R0tV = R0'*V;
alpha = alpha - (RV'*(RV*alpha)).*dd;     % long alpha vector for ordinary infEP
post.alpha = R0tV*alpha;       % alpha = R0'*V*inv(Kt+diag(1./ttau))*(tnu./ttau)
B = R0tV.*repmat(dd',nu,1); L = B*R0tV'; B = B*RV';
post.L = B*B' - L;                      % L = -R0'*V*inv(Kt+diag(1./ttau))*V'*R0

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  RVdd = RV.*repmat(dd',nu,1);
  for i=1:length(hyp.cov)
    [ddiagK,dKuu,dKu] = feval(cov{:}, hyp.cov, x, [], i); % eval cov derivatives
    dA = 2*dKu'-R0tV'*dKuu;                                       % dQ = dA*R0tV
    w = sum(dA.*R0tV',2); v = ddiagK-w;   % w = diag(dQ); v = diag(dK)-diag(dQ);
    z = dd'*(v+w) - sum(RVdd.*RVdd,1)*v - sum(sum( (RVdd*dA)'.*(R0tV*RVdd') ));
    dnlZ.cov(i) = (z - alpha'*(alpha.*v) - (alpha'*dA)*(R0tV*alpha))/2;
  end
  for i = 1:numel(hyp.lik)                                   % likelihood hypers
    dlik = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf, i);
    dnlZ.lik(i) = -sum(dlik);
    if i==numel(hyp.lik)
      % since snu2 is a fixed fraction of sn2, there is a covariance-like term
      % in the derivative as well
      v = sum(R0tV.*R0tV,1)';
      z = sum(sum( (RVdd*R0tV').^2 )) - sum(RVdd.*RVdd,1)*v;
      z = z + post.alpha'*post.alpha - alpha'*(v.*alpha);
      dnlZ.lik(i) = dnlZ.lik(i) + snu2*z;
    end
  end
  [junk,dlZ] = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);% mean hyps
  for i = 1:numel(hyp.mean)
    dm = feval(mean{:}, hyp.mean, x, i);
    dnlZ.mean(i) = -dlZ'*dm;
  end
  if isfield(hyp,'xu')                   % derivatives w.r.t. inducing points xu
    xu = cov{3};
    cov = cov{2};             % get the non FITC part of the covariance function
    Kpu  = cov_deriv_sq_dist(cov,hyp.cov,xu,x);             % d K(xu,x ) / d D^2
    Kpuu = cov_deriv_sq_dist(cov,hyp.cov,xu);               % d K(xu,xu) / d D^2
    if iscell(cov), covstr = cov{1}; else covstr = cov; end
    if ~ischar(covstr), covstr = func2str(covstr); end
    if numel(strfind(covstr,'iso'))>0              % characteristic length scale
      e = 2*exp(-2*hyp.cov(1));
    else
      e = 2*exp(-2*hyp.cov(1:D));
    end
    B = (R0'*R0)*Ku;

    W = ttau;
    t = W./(1+W.*d0);
    diag_dK = alpha.*alpha + sum(RVdd.*RVdd,1)' - t;
    v = diag_dK+t;                 % BdK = B * ( dnlZ/dK - diag(diag(dnlZ/dK)) )
    BdK = (B*alpha)*alpha' - B.*repmat(v',nu,1);
    BdK = BdK + (B*RVdd')*RVdd;
    A = Kpu.*BdK; C = Kpuu.*(BdK*B'); C = diag(sum(C,2)-sum(A,2)) - C;
    dnlZ.xu = A*x*diag(e) + C*xu*diag(e);    % bring in data and inducing points
  end
end

% refresh the representation of the posterior from initial and site parameters
% to prevent possible loss of numerical precision after many epfitcUpdates
% effort is O(n*nu^2) provided that nu<n
function [d,P,R,nn,gg] = epfitcRefresh(d0,P0,R0,R0P0, w,b)
  nu = size(R0,1);                                   % number of inducing points
  rot180   = @(A)   rot90(rot90(A));                   % little helper functions
  chol_inv = @(A) rot180(chol(rot180(A))')\eye(nu);               % chol(inv(A))
  t  = 1./(1+d0.*w);                                   % temporary variable O(n)
  d  = d0.*t;                                                             % O(n)
  P  = repmat(t',nu,1).*P0;                                            % O(n*nu)
  T  = repmat((w.*t)',nu,1).*R0P0;                % temporary variable O(n*nu^2)
  R  = chol_inv(eye(nu)+R0P0*T');                                    % O(n*nu^3)
  nn = d.*b;                                                              % O(n)
  gg = R0'*(R'*(R*(R0P0*(t.*b))));                                     % O(n*nu)

% compute the marginal likelihood approximation
% effort is O(n*nu^2) provided that nu<n
function [nlZ,nu_n,tau_n] = ...
                        epfitcZ(d,P,R,nn,gg,ttau,tnu, d0,R0,P0, y,lik,hyp,m,inf)
  T = (R*R0)*P;                                             % temporary variable
  diag_sigma = d + sum(T.*T,1)'; mu = nn + P'*gg;       % post moments O(n*nu^2)
  tau_n = 1./diag_sigma-ttau;              % compute the log marginal likelihood
  nu_n  = mu./diag_sigma-tnu+m.*tau_n;            % vectors of cavity parameters
  lZ = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, inf);
  nu = size(gg,1);
  U = (R0*P0)'.*repmat(1./sqrt(d0+1./ttau),1,nu);
  L = chol(eye(nu)+U'*U);
  ld = 2*sum(log(diag(L))) + sum(log(1+d0.*ttau));
  t = T*tnu; tnu_Sigma_tnu = tnu'*(d.*tnu) + t'*t;
  nlZ = ld/2 -sum(lZ) -tnu_Sigma_tnu/2  ...
    -(nu_n-m.*tau_n)'*((ttau./tau_n.*(nu_n-m.*tau_n)-2*tnu)./(ttau+tau_n))/2 ...
      +sum(tnu.^2./(tau_n+ttau))/2-sum(log(1+ttau./tau_n))/2;

% update the representation of the posterior to reflect modification of the site 
% parameters by w(i) <- wi and b(i) <- bi
% effort is O(nu^2)
% Pi = P(:,i) is passed instead of P to prevent allocation of a new array
function [d,Pi,R,nn,gg,w,b] = ...
                           epfitcUpdate(d,Pi,R,nn,gg, w,b, i,wi,bi, m,d0,P0,R0)
  dwi = wi-w(i); dbi = bi-b(i);
  hi = nn(i) + m(i) + Pi'*gg;                   % posterior mean of site i O(nu)
  t = 1+dwi*d(i);
  d(i) = d(i)/t;                                                          % O(1)
  nn(i) = d(i)*bi;                                                        % O(1)
  r = 1+d0(i)*w(i);
  r = r*r/dwi + r*d0(i);
  v = R*(R0*P0(:,i));
  r = 1/(r+v'*v);
  if r>0
    R = cholupdate(R,sqrt( r)*R'*v,'-');
  else
    R = cholupdate(R,sqrt(-r)*R'*v,'+');
  end
  gg = gg + ((dbi-dwi*(hi-m(i)))/t)*(R0'*(R'*(R*(R0*Pi))));            % O(nu^2)
  w(i) = wi; b(i) = bi;                            % update site parameters O(1)
  Pi = Pi/t;                                                             % O(nu)
