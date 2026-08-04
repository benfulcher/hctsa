function [post nlZ dnlZ] = infState(hyp, mean, cov, lik, x, y, opt)

% Inference for a GP with state space covariance.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters.
% See also "help infMethods".
%
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with apxState.m.
%
% Copyright (c) by Hannes Nickisch 2018-05-30.
%
% See also infMethods.m, cov/apxState.m.

if nargin<7, opt = []; end                          % make sure parameter exists

if ischar(lik) || isa(lik,'function_handle'), lik = {lik};  end      % make cell
if iscell(lik), lstr = lik{1}; else lstr = lik; end
if isa(lstr,'function_handle'), lstr = func2str(lstr); end
if ischar(cov) || isa(cov,'function_handle'), cov  = {cov};  end     % make cell
cstr = cov{1}; if isa(cstr,'function_handle'), cstr = func2str(cstr); end
if isequal(cstr,'apxState'), cov = cov{2}; end
if nargout==0, return, end                                  % nothing to do here
if ~isfield(opt,'inf')                             % default is GaussLik/Laplace
  if isequal(lstr,'likGauss'), inf = @infGaussLik; else inf = @infLaplace; end
else inf = opt.inf;
end
if ischar(inf) || isa(inf,'function_handle'), inf = {inf};  end      % make cell
if iscell(inf), istr = inf{1}; else istr = inf; end
if isa(istr,'function_handle'), istr = func2str(istr); end

[m,dm] = feval(mean{:},hyp.mean,x);                     % evaluate mean function
if strcmp('infEP',istr)
  mom = @(mu,s2,l) feval(lik{:},hyp.lik,y(l),mu+m(l),s2,'infEP');      % for ADF
  n = size(x,1); tnu = zeros(n,1); ttau = zeros(n,1);
  [fmu,fs2,ymu,ys2,z,tnu,ttau] = predict(x,hyp,mean,cov,lik,x,tnu,ttau,mom);
  post.sW = sqrt(ttau);
  post.tnu = tnu;
  K = apx(hyp,{@apxState,cov},x,opt);
  [ldB2,solveKiW,dW,dhyp] = K.fun(ttau);              % functions depending on W
  diagSigma = 2*dW;
  post.alpha = tnu-solveKiW(K.mvm(tnu)+m);
  mu = K.mvm(post.alpha)+m;
  tau_n = 1./diagSigma-ttau; nu_n = mu./diagSigma-tnu;                % cavities
  post.L = @(r)-K.P(solveKiW(K.Pt(r)));
  lZ = feval(lik{:}, hyp.lik, y, nu_n./tau_n, 1./tau_n, 'infEP');
  p = tnu-m.*ttau; q = nu_n-m.*tau_n; r = K.mvm(p);          % auxiliary vectors
  nlZ = ldB2-sum(lZ)-p'*r/2 +r'*solveKiW(r)/2 +(diagSigma'*p.^2)/2 ...
              -q'*((ttau./tau_n.*q-2*p).*diagSigma)/2-sum(log(1+ttau./tau_n))/2;
  if nargout>2
    dnlZ = dhyp(post.alpha);                         % covariance-related hypers
    dnlZ.lik = zeros(1,numel(hyp.lik));                        % allocate memory
    for i = 1:numel(hyp.lik)                                 % likelihood hypers
      dlik = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,inf,i);
      dnlZ.lik(i) = -sum(dlik);
    end
    [junk,dlZ] = feval(lik{:},hyp.lik,y,nu_n./tau_n,1./tau_n,'infEP');    % mean
    dnlZ.mean = -dm(dlZ);
  end
else
  if nargout<3
    [post nlZ]      = feval(inf{:}, hyp, mean, {@apxState,cov}, lik, x, y, opt);
  else
    [post nlZ dnlZ] = feval(inf{:}, hyp, mean, {@apxState,cov}, lik, x, y, opt);
  end
  W = post.sW.^2;
  if isequal(lstr,'likGauss')                               % compute 'residual'
    post.tnu = W.*(y-m);                         % all right for pure regression
  else
    K = apx(hyp,{@apxState,cov},x,opt);
    post.tnu = W.*K.mvm(post.alpha) + post.alpha;
  end
end

% add fast predictions to post structure, f|y,mu|s2
post.predict = @(xs) predict(xs,hyp,mean,cov,lik,x,post.tnu,post.sW.^2,[]);

% Compute latent and predictive means and variances by Kalman filtering.
function [fmu,fs2,ymu,ys2,z,tnu,ttau] = ...
                                     predict(xs,hyp,mean,cov,lik,x,tnu,ttau,mom)
  [F,L,Qc,H,Pinf,stat] = apxState(cov,hyp.cov);       % set up state space model
  [As,Qs, t,tid] = apxState('trans',F,L,Qc,H,Pinf,stat,[x;xs]);
  n = size(x,1); ns = size(xs,1);   % t(id(1:n))==x, t(id(n+(1:ns)))==t(sid)==xs
  % xx = zeros(n+ns+1,1); xx(tid) = t; => xx(1:n)==x, xx(n+(1:ns))==xs
  id = (1:n+ns+1)'; id(tid) = id; sid = id(n+(1:ns));
  [z,ga,tnu,ttau,ms,Ps] = apxState('kalman',As,Qs, H,Pinf, tnu,ttau, t,tid,mom);
  [ms,Ps] = apxState('rts',As,Qs,ms,Ps);         % Rauch-Tung-Striebel smoothing
  fmu = (H*ms(:,sid))' + feval(mean{:},hyp.mean,xs);
  fs2 = arrayfun(@(k) H*Ps(:,:,k)*H', sid);
  [lp,ymu,ys2] = feval(lik{:},hyp.lik,[],fmu,fs2);% propagate through likelihood