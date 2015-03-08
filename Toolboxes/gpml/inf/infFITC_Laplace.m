function [post nlZ dnlZ] = infFITC_Laplace(hyp, mean, cov, lik, x, y)

% FITC-Laplace approximation to the posterior Gaussian process. The function is
% equivalent to infLaplace with the covariance function:
%   Kt = Q + G; G = diag(g); g = diag(K-Q);  Q = Ku'*inv(Kuu + snu2*eye(nu))*Ku;
% where Ku and Kuu are covariances w.r.t. to inducing inputs xu and
% snu2 = sn2/1e6 is the noise of the inducing inputs. We fixed the standard
% deviation of the inducing inputs snu to be a one per mil of the measurement 
% noise's standard deviation sn. In case of a likelihood without noise
% parameter sn2, we simply use snu2 = 1e-6.
%
% The implementation exploits the Woodbury matrix identity
%   inv(Kt) = inv(G) - inv(G)*Ku'*inv(Kuu+Ku*inv(G)*Ku')*Ku*inv(G)
% in order to be applicable to large datasets. The computational complexity
% is O(n nu^2) where n is the number of data points x and nu the number of
% inducing inputs in xu.
% The posterior N(f|h,Sigma) is given by h = m+mu with mu = nn + P'*gg and
% Sigma = inv(inv(K)+diag(W)) = diag(d) + P'*R0'*R'*R*R0*P.
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
% Copyright (c) by Hannes Nickisch, 2013-10-28.
%
% See also INFMETHODS.M, COVFITC.M.

persistent last_alpha                                   % copy of the last alpha
if any(isnan(last_alpha)), last_alpha = zeros(size(last_alpha)); end   % prevent
tol = 1e-6;                   % tolerance for when to stop the Newton iterations
smax = 2; Nline = 10; thr = 1e-4;                       % line search parameters
maxit = 20;                                    % max number of Newton steps in f
inf = 'infLaplace';
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

Psi_old = Inf;    % make sure while loop starts by the largest old objective val
if any(size(last_alpha)~=[n,1])     % find a good starting point for alpha and f
  alpha = zeros(n,1); f = mvmK(alpha,V,d0)+m; % start at mean if sizes not match                  
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp; Psi_new = -sum(lp);
else
  alpha = last_alpha; f = mvmK(alpha,V,d0)+m;                     % try last one
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp;
  Psi_new = alpha'*(f-m)/2 - sum(lp);                 % objective for last alpha
  Psi_def = -feval(lik{:},hyp.lik,y,m,[],inf); % objective for default init f==m
  if Psi_def < Psi_new                         % if default is better, we use it
    alpha = zeros(n,1); f = mvmK(alpha,V,d0)+m;
    [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp; Psi_new=-sum(lp);
  end
end
isWneg = any(W<0);       % flag indicating whether we found negative values of W
it = 0;                            % this happens for the Student's t likelihood

while Psi_old - Psi_new > tol && it<maxit                         % begin Newton
  Psi_old = Psi_new; it = it+1;
  if isWneg       % stabilise the Newton direction in case W has negative values
    W = max(W,0);      % stabilise the Hessian to guarantee postive definiteness
    tol = 1e-8;            % increase accuracy to also get the derivatives right
    % In Vanhatalo et. al., GPR with Student's t likelihood, NIPS 2009, they use
    % a more conservative strategy then we do being equivalent to 2 lines below.
    % nu  = exp(hyp.lik(1));                  % degree of freedom hyperparameter
    % W  = W + 2/(nu+1)*dlp.^2;               % add ridge according to Vanhatalo
  end
  b = W.*(f-m) + dlp; dd = 1./(1+W.*d0);
  RV = chol_inv(eye(nu)+(V.*repmat((W.*dd)',nu,1))*V')*V;
  dalpha = dd.*b - (W.*dd).*(RV'*(RV*(dd.*b))) - alpha; % Newt dir + line search
  [s,Psi_new,Nfun,alpha,f,dlp,W] = brentmin(0,smax,Nline,thr, ...
                                 @Psi_line,4,dalpha,alpha,hyp,V,d0,m,lik,y,inf);
  isWneg = any(W<0);
end                                                    % end Newton's iterations

last_alpha = alpha;                                     % remember for next call
[lp,dlp,d2lp,d3lp] = feval(lik{:},hyp.lik,y,f,[],inf); W=-d2lp; isWneg=any(W<0);
post.alpha = R0'*(V*alpha);                    % return the posterior parameters
post.sW = sqrt(abs(W)).*sign(W);             % preserve sign in case of negative
dd = 1./(1+d0.*W);                                     % temporary variable O(n)
A = eye(nu)+(V.*repmat((W.*dd)',nu,1))*V';        % temporary variable O(n*nu^2)
R0tV = R0'*V; B = R0tV.*repmat((W.*dd)',nu,1);   % temporary variables O(n*nu^2)
post.L = -B*R0tV';          % L = -R0'*V*inv(Kt+diag(1./ttau))*V'*R0, first part
if any(1+d0.*W<0)
  B = B*V'; post.L = post.L + (B*inv(A))*B';
  nlZ = NaN; dnlZ = struct('cov',0*hyp.cov, 'mean',0*hyp.mean, 'lik',0*hyp.lik);
  warning('W is too negative; nlZ and dnlZ cannot be computed.')
  return
end
nlZ = alpha'*(f-m)/2 - sum(lp) - sum(log(dd))/2 + sum(log(diag(chol(A))));
RV = chol_inv(A)*V; RVdd = RV.*repmat((W.*dd)',nu,1);     % RVdd needed for dnlZ
B = B*RV'; post.L = post.L + B*B';

if nargout>2                                           % do we want derivatives?
  dnlZ = hyp;                                   % allocate space for derivatives
  [d,P,R] = fitcRefresh(d0,Ku,R0,V,W);               % g = diag(inv(inv(K)+W))/2
  g = d/2 + sum(((R*R0)*P).^2,1)'/2;
  t = W./(1+W.*d0);

  dfhat = g.*d3lp;  % deriv. of nlZ wrt. fhat: dfhat=diag(inv(inv(K)+W)).*d3lp/2
  for i=1:length(hyp.cov)                                    % covariance hypers
    [ddiagK,dKuu,dKu] = feval(cov{:}, hyp.cov, x, [], i); % eval cov derivatives
    dA = 2*dKu'-R0tV'*dKuu;                                       % dQ = dA*R0tV
    w = sum(dA.*R0tV',2); v = ddiagK-w;   % w = diag(dQ); v = diag(dK)-diag(dQ);
    dnlZ.cov(i) = ddiagK'*t - sum(RVdd.*RVdd,1)*v;               % explicit part
    dnlZ.cov(i) = dnlZ.cov(i)-sum(sum((RVdd*dA).*(RVdd*R0tV'))); % explicit part
    dnlZ.cov(i) = dnlZ.cov(i)/2-alpha'*(dA*(R0tV*alpha)+v.*alpha)/2;  % explicit
    b = dA*(R0tV*dlp) + v.*dlp;            % b-K*(Z*b) = inv(eye(n)+K*diag(W))*b
    KZb = mvmK(mvmZ(b,RVdd,t),V,d0);
    dnlZ.cov(i) = dnlZ.cov(i) - dfhat'*( b-KZb );                % implicit part
  end
  for i=1:length(hyp.lik)                                    % likelihood hypers
    [lp_dhyp,dlp_dhyp,d2lp_dhyp] = feval(lik{:},hyp.lik,y,f,[],inf,i);
    dnlZ.lik(i) = -g'*d2lp_dhyp - sum(lp_dhyp);                  % explicit part
    b = mvmK(dlp_dhyp,V,d0);                                     % implicit part
    dnlZ.lik(i) = dnlZ.lik(i) - dfhat'*(b-mvmK(mvmZ(b,RVdd,t),V,d0));
    if i==numel(hyp.lik)
      % since snu2 is a fixed fraction of sn2, there is a covariance-like term
      % in the derivative as well
      snu = sqrt(snu2);
      T = chol_inv(Kuu + snu2*eye(nu)); T = T'*(T*(snu*Ku)); t = sum(T.*T,1)';
      z = alpha'*(T'*(T*alpha)-t.*alpha) - sum(RVdd.*RVdd,1)*t;
      z = z + sum(sum( (RVdd*T').^2 ));
      b = (t.*dlp-T'*(T*dlp))/2;
      KZb = mvmK(mvmZ(b,RVdd,t),V,d0);
      z = z - dfhat'*( b-KZb );
      dnlZ.lik(i) = dnlZ.lik(i) + z;
    end
  end
  for i=1:length(hyp.mean)                                         % mean hypers
    dm = feval(mean{:}, hyp.mean, x, i);
    dnlZ.mean(i) = -alpha'*dm;                                   % explicit part
    Zdm = mvmZ(dm,RVdd,t);
    dnlZ.mean(i) = dnlZ.mean(i) - dfhat'*(dm-mvmK(Zdm,V,d0));    % implicit part
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
    q = dfhat - mvmZ(mvmK(dfhat,V,d0),RVdd,t);                   % implicit part
    B = (R0'*R0)*Ku;
    diag_dK = alpha.*alpha + sum(RVdd.*RVdd,1)' - t + 2*dlp.*q;
    v = diag_dK+t;                 % BdK = B * ( dnlZ/dK - diag(diag(dnlZ/dK)) )
    BdK = (B*alpha)*alpha' - B.*repmat(v',nu,1) + (B*dlp)*q' + (B*q)*dlp';
    BdK = BdK + (B*RVdd')*RVdd;
    A = Kpu.*BdK; C = Kpuu.*(BdK*B'); C = diag(sum(C,2)-sum(A,2)) - C;
    dnlZ.xu = A*x*diag(e) + C*xu*diag(e);    % bring in data and inducing points
  end
end

% matrix vector multiplication with Z=inv(K+inv(W))
function Zx = mvmZ(x,RVdd,t)
  Zx = t.*x - RVdd'*(RVdd*x);

% matrix vector multiplication with approximate covariance matrix
function Kal = mvmK(al,V,d0)
  Kal = V'*(V*al) + d0.*al;

% criterion Psi at alpha + s*dalpha for line search
function [Psi,alpha,f,dlp,W] = Psi_line(s,dalpha,alpha,hyp,V,d0,m,lik,y,inf)
  alpha = alpha + s*dalpha;
  f = mvmK(alpha,V,d0)+m;
  [lp,dlp,d2lp] = feval(lik{:},hyp.lik,y,f,[],inf); W = -d2lp;
  Psi = alpha'*(f-m)/2 - sum(lp);

% refresh the representation of the posterior from initial and site parameters
% to prevent possible loss of numerical precision after many epfitcUpdates
% effort is O(n*nu^2) provided that nu<n
% Sigma = inv(inv(K)+diag(W)) = diag(d) + P'*R0'*R'*R*R0*P.
function [d,P,R] = fitcRefresh(d0,P0,R0,R0P0, w)
  nu = size(R0,1);                                   % number of inducing points
  rot180   = @(A)   rot90(rot90(A));                   % little helper functions
  chol_inv = @(A) rot180(chol(rot180(A))')\eye(nu);               % chol(inv(A))
  t  = 1./(1+d0.*w);                                   % temporary variable O(n)
  d  = d0.*t;                                                             % O(n)
  P  = repmat(t',nu,1).*P0;                                            % O(n*nu)
  T  = repmat((w.*t)',nu,1).*R0P0;                % temporary variable O(n*nu^2)
  R  = chol_inv(eye(nu)+R0P0*T');                                    % O(n*nu^3)
