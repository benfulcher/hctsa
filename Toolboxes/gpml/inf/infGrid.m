function [post nlZ dnlZ] = infGrid(hyp, mean, cov, lik, x, y, opt)

% Inference for a GP with grid-based approximate covariance.
%
% The (Kronecker) covariance matrix used is given by:
%   K = kron( kron(...,K{2}), K{1} ) = K_p x .. x K_2 x K_1.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters.
% The result is exact for complete grids, otherwise results are approximate.
% See also "help infMethods".
%
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with apxGrid.m.
%
% In case of equispaced data points, we use Toeplitz/BTTB algebra. We use a
% circulant embedding approach to approximate the log determinant of the
% covariance matrix. If any of the factors K_i, i=1..p has Toeplitz or more
% general BTTB structure (which is indicated by K.kron.factor(i).descr being
% equal to 'toep', 'bttb2', 'bttb3', etc.), we automatically use the circulant
% determinant approximation. The grid specification needs to reflect this.
% There are some examples to illustrate the doubly nested curly bracket
% formalism. See also "help apxGrid".
%
% There are a set of options available:
% opt.pred_var, minimum value is 20 as suggested in the Papandreou paper
%   Instead of the data x, we can tell the engine to use x*hyp.P' to make grid
%   methods available to higher dimensional data. We offer two ways of
%   restricting the projection matrix hyp.P to either orthonormal matrices,
%   where hyp.P*hyp.P'=I or normalised projections diag(hyp.P*hyp.P')=1.
%   Negative values indicate that we are using Lanczos variance estimation
%   as described by Pleiss et al. in https://arxiv.org/abs/1803.06058.
% opt.proj = 'orth'; enforce orthonormal projections by working with 
%   sqrtm(hyp.P*hyp.P')\hyp.P instead of hyp.P
% opt.proj = 'norm'; enforce normal projections by working with 
%   diag(1./sqrt(diag(hyp.P*hyp.P')))*hyp.P instead of hyp.P
%
% There are a number of options inherited from apx.m
%   The conjugate gradient-based linear system solver has two adjustable
%   parameters, the relative residual threshold for convergence opt.cg_tol and
%   the maximum number of MVMs opt.cg_maxit until the process stops.
%  opt.cg_tol,   default is 1e-6      as in Matlab's pcg function
%  opt.cg_maxit, default is min(n,20) as in Matlab's pcg function
%    We can tell the inference engine to make functions post.fs2 and post.ys2
%    available in order to compute the latent and predictive variance of an
%    unknown test data point. Precomputations for Perturb-and-MAP sampling are
%    required for these functions.
%  opt.stat = true returns a little bit of output summarising the exploited
%    structure of the covariance of the grid.
%    Please see cov/apxGrid.m for details.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2018-03-26.
%
% See also infMethods.m, cov/apxGrid.m.

if nargin<7, opt = []; end                          % make sure parameter exists
xg = cov{3}; p = numel(xg);                 % extract underlying grid parameters
[ng,Dg] = apxGrid('size',xg); N = prod(ng); D = sum(Dg);        % dimensionality
if isfield(opt,'proj'), proj = opt.proj; else proj = 'none'; end    % projection
hP = 1;                                                   % no projection at all
if isfield(hyp,'proj')                 % apply transformation matrix if provided
  hP = hyp.proj;
  if strncmpi(proj,'orth',4)
    hP = sqrtm(hP*hP'+eps*eye(D))\hP;                    % orthonormal projector
  elseif strncmpi(proj,'norm',4)
    hP = diag(1./sqrt(diag(hP*hP')+eps))*hP;                  % normal projector
  end
end
if isfield(opt,'deg'), deg = opt.deg; else deg = 3; end       % interpol. degree
% no of samples for covariance hyperparameter sampling approximation,
% see Po-Ru Loh et.al.: "Contrasting regional architectures of schizophrenia and
% other complex diseases using fast variance components analysis, biorxiv.org
if isfield(opt,'ndcovs'), ndcovs = max(opt.ndcovs,20);
else ndcovs = 0; end
[K,M] = feval(cov{:}, hyp.cov, x*hP');    % evaluate covariance mat constituents
m = feval(mean{:}, hyp.mean, x*hP');                      % evaluate mean vector
if iscell(lik), lstr = lik{1}; else lstr = lik; end
if isa(lstr,'function_handle'), lstr = func2str(lstr); end
if isequal(lstr,'likGauss'), inf = @infGaussLik; else inf = @infLaplace; end
if nargout>0
  if nargout<3
    [post nlZ] = inf(hyp, mean, cov, lik, x*hP', y, opt);
  else
    [post nlZ dnlZ] = inf(hyp, mean, cov, lik, x*hP', y, opt);
    if isfield(hyp,'proj')
      dnlZ.proj=deriv_proj(post.alpha,hP,K,covGrid('flatten',xg),m,mean,hyp,x);
      if     strncmpi(proj,'orth',4), dnlZ.proj=chain_orth(hyp.proj,dnlZ.proj);
      elseif strncmpi(proj,'norm',4), dnlZ.proj=chain_norm(hyp.proj,dnlZ.proj);
      end
    end
  end
else return, end

% no of samples for perturb-and-MAP, see George Papandreou and Alan L. Yuille:
% "Efficient Variational Inference in Large-Scale Bayesian Compressed Sensing"
ns = 0;                       % do nothing per default, 20 is suggested in paper
if isfield(opt,'pred_var')                              % at least default value
  ns = sign(opt.pred_var) * max(ceil(abs(opt.pred_var)),20);
end
if ndcovs>0 && nargout>2                            % possibly draw more samples
  ns = sign(ns) * max(abs(ns),ndcovs);
end  
Mtal = M'*post.alpha;                         % blow up alpha vector from n to N
kronmvm = K.kronmvm;
if ns>0
  s = 3;                                      % Whittle embedding overlap factor
  [V,ee,e] = apxGrid('eigkron',K,xg,s);            % perform eigen-decomposition
  % explained variance on the grid vg=diag(Ku*M'*inv(C)*M*Ku), C=M*Ku*M'+inv(W)
  % relative accuracy r = std(vg_est)/vg_exact = sqrt(2/ns)
  A = sample(V,e,M,post.sW,ns,kronmvm); A = post.L(A);           % a~N(0,inv(C))
  z = K.mvm(M'*A); vg = sum(z.*z,2)/ns;             % z ~ N(0,Ku*M'*inv(C)*M*Ku)
  if ndcovs>0
    dnlZ.covs = - apxGrid('dirder',K,xg,M,post.alpha,post.alpha)/2;
    na = size(A,2);
    for i=1:na                                % compute (E[a'*dK*a] - a'*dK*a)/2
      dnlZ.covs = dnlZ.covs + apxGrid('dirder',K,xg,M,A(:,i),A(:,i))/(2*na);
    end
    if isfield(hyp,'proj')                              % data projection matrix
      dPs = zeros(size(hyp.proj));                  % allocate memory for result
      KMtal = K.mvm(Mtal); [M,dM]=covGrid('interp',xg,x*hP');       % precompute
      for i = 1:size(dPs,1)
        if equi(xg,i), wi = max(xg{i})-min(xg{i}); else wi = 1; end    % scaling
        for j = 1:size(dPs,2)
            dMtal = dM{i}'*(x(:,j).*post.alpha/wi);
            dMtA  = dM{i}'*(repmat(x(:,j),1,na).*A/wi);
            dPs(i,j) = sum(sum(dMtA.*z))/ndcovs - dMtal'*KMtal;
        end
      end
      if     strncmpi(proj,'orth',4),dnlZ.projs = chain_orth(hyp.proj,dPs);
      elseif strncmpi(proj,'norm',4),dnlZ.projs = chain_norm(hyp.proj,dPs);
      else                           dnlZ.projs = dPs; end
    end
  end
elseif ns<0                                       % variance estimate using LOVE
  v = M*K.mvm( ones(N,1)/N ); n = size(M,1); a = 1./post.sW.^2;
  [Q,T] = apxGrid('lanczos_arpack', @(x)M*K.mvm(M'*x)+a.*x, v,min(abs(ns),n-1));
  R = K.mvm(M'*Q)/chol(T); vg = sum(R.^2,2);
else
  vg = zeros(N,1);                                       % no variance explained
end
% add fast predictions to post structure, f|y,mu|s2
post.predict = @(xs) predict(xs*hP',xg,K.mvm(Mtal),vg,hyp,mean,cov,lik,deg);
% global mem, S = whos(); mem=0; for i=1:numel(S), mem=mem+S(i).bytes/1e6; end

% Compute latent and predictive means and variances by grid interpolation.
function [fmu,fs2,ymu,ys2] = predict(xs,xg,Kalpha,vg,hyp,mean,cov,lik,deg)
  Ms = apxGrid('interp',xg,xs,deg);                % obtain interpolation matrix
  xs = apxGrid('idx2dat',xg,xs,deg);                    % deal with index vector
  ms = feval(mean{:},hyp.mean,xs);                         % evaluate prior mean
  fmu = ms + Ms*Kalpha;                 % combine and perform grid interpolation
  if nargout>1
    if norm(vg,1)>1e-10, ve = Ms*vg; else ve = 0; end    % interp grid var expl.
    ks = feval(cov{:},hyp.cov,xs,'diag');              % evaluate prior variance
    fs2 = max(ks-ve,0);              % combine, perform grid interpolation, clip
    if nargout>2, [lp, ymu, ys2] = feval(lik{:},hyp.lik,[],fmu,fs2); end
  end

% sample a~N(0,C), C = M*Ku*M'+inv(W), W=sW^2
function A = sample(V,e,M,sW,ns,kronmvm)
  [n,N] = size(M);
  A = randn(N,ns);                                                    % a~N(0,I)
  A = kronmvm(V,repmat(sqrt(e),1,ns).*kronmvm(V,A,1));               % a~N(0,Ku)
  A = M*A + bsxfun(@times,1./sW,randn(n,ns));

% compute derivative of neg log marginal likelihood w.r.t. projection matrix P
function dP = deriv_proj(alpha,P,K,xg,m,mean,hyp,x)
  xP = x*P'; [M,dM] = covGrid('interp',xg,xP); % grid interp derivative matrices
  beta = K.mvm(M'*alpha);                          % dP(i,j) = -alpha'*dMij*beta
  dP = zeros(size(P)); h = 1e-4;               % allocate result, num deriv step
  for i=1:size(P,1)
    if equi(xg,i), wi = max(xg{i})-min(xg{i}); else wi = 1; end % scaling factor
    xP(:,i) = xP(:,i)+h;            % perturb ith component of projection matrix
    dmi = (feval(mean{:},hyp.mean,xP)-m)/h;    % numerically estimate dm/dP(:,i)
    xP(:,i) = xP(:,i)-h;                                     % undo perturbation
    betai = dmi + dM{i}*beta/wi;
    for j=1:size(P,2), dP(i,j) = -alpha'*(x(:,j).*betai); end
  end

function eq = equi(xg,i)                        % grid along dim i is equispaced
  ni = size(xg{i},1);
  if ni>1                              % diagnose if data is linearly increasing
    dev = abs(diff(xg{i})-ones(ni-1,1)*(xg{i}(2,:)-xg{i}(1,:)));
    eq = max(dev(:))<1e-9;
  else
    eq = true;
  end

% chain rule for the function Q = sqrtm(P*P')\P;  for d sqrtm(X) see the website
function dQ = chain_orth(P,dP)  % http://math.stackexchange.com/questions/540361
  [V,F] = eig(P*P'); sf = sqrt(diag(F)); S = V*diag(sf)*V';         % eig-decomp
  H = dP'/S; G = H'*(P'/S); o = ones(size(dP,1),1);                 % chain rule
  dQ = (H - P'*V*((V'*(G+G')*V)./(sf*o'+o*sf'))*V')';

% chain rule for the function Q = diag(1./sqrt(diag(P*P')))*P;
function dQ = chain_norm(P,dP)
  p = 1./sqrt(diag(P*P'));
  dQ = diag(p)*dP - diag(diag(dP*P').*p.^3)*P;