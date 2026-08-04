function [post, nlZ, dnlZ] = infVB(hyp, mean, cov, lik, x, y, opt)

% Variational approximation to the posterior Gaussian process.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Minimisation of an upper bound on the negative marginal likelihood using a
% sequence of infLaplace calls where the smoothed likelihood
% likVB(f) = lik(..,g,..) * exp(b*(f-g)), g = sign(f-z)*sqrt((f-z)^2+v)+z, where
%     v   .. marginal variance = (positive) smoothing width, and
%     lik .. lik function such that p(y|f)=lik(..,f,..).
%
% The problem is convex whenever the likelihood is log-concave. At the end, the
% optimal width W is obtained analytically.
%
% Copyright (c) by Hannes Nickisch 2016-10-21.
%
% See also infMethods.m, gp.m.

n = size(x,1);
if nargin<=6, opt = []; end                        % make opt variable available
if isstruct(cov), K = cov;                   % use provided covariance structure
else K = apx(hyp,cov,x,opt); end               % set up covariance approximation
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv
if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end

if nargin<=6, opt = []; end                        % make opt variable available
if ~isfield(opt,'postL'), opt.postL = false; end   % not compute L in infLaplace
if isfield(opt,'out_nmax'), out_nmax = opt.out_nmax; % maximal no of outer loops
else out_nmax = 15; end                                          % default value
if isfield(opt,'out_tol'),  out_tol  = opt.out_tol;     % outer loop convergence
else out_tol = 1e-5; end                                         % default value

sW = ones(n,1);                                % init with some reasonable value
opt.postL = false;                    % avoid computation of L inside infLaplace
for i=1:out_nmax
  [junk,junk,dW] = K.fun(sW.*sW); v = 2*dW;
  [post,junk,junk,alpha] = infLaplace(hyp, mean, K, {@likVB,v,lik}, x, y, opt);
  % post.sW is very different from the optimal sW for non Gaussian likelihoods
  sW_old = sW; f = K.mvm(alpha)+m;                              % posterior mean
  [lp,junk,junk,sW,b,z] = feval(@likVB,v,lik,hyp.lik,y,f);
  if max(abs(sW-sW_old))<out_tol, break, end              % diagnose convergence
end

post.sW = sW;                                             % posterior parameters
[ldB2,solveKiW,dW,dhyp] = K.fun(sW.*sW); post.L = @(r) -K.P(solveKiW(K.Pt(r)));

ga = 1./(sW.*sW); be = b+z./ga;        % variance, lower bound offset from likVB
h = f.*(2*be-f./ga) - 2*lp - v./ga;     % h(ga) = s*(2*b-f/ga)+ h*(s) - v*(1/ga)
t = b.*ga+z-m; nlZ = ldB2 + (sum(h)+t'*solveKiW(t)-(be.*be)'*ga )/2; % var bound

if nargout>2                                           % do we want derivatives?
  dnlZ = dhyp(alpha);                                        % covariance hypers
  dnlZ.lik = zeros(size(hyp.lik));                             % allocate memory
  if ~strcmp(likstr,'likGauss')                              % likelihood hypers
    for j=1:length(hyp.lik)
      sign_fmz = 2*(f-z>=0)-1;                % strict sign mapping; sign(0) = 1
      g = sign_fmz.*sqrt((f-z).^2 + v) + z;
      dhhyp = -2*feval(lik{:},hyp.lik,y,g,[],'infLaplace',j);
      dnlZ.lik(j) = sum(dhhyp)/2;
    end  
  else                                 % special treatment for the Gaussian case
    sn2 = exp(2*hyp.lik); dnlZ.lik = -sn2*(alpha'*alpha) - 2*sum(dW)/sn2 + n;
  end
  dnlZ.mean = -dm(alpha);                                          % mean hypers
end

% Smoothed likelihood function; instead of p(y|f)=lik(..,f,..) compute
%   likVB(f) = lik(..,g,..)*exp(b*(f-g)), g = sign(f-z)*sqrt((f-z)^2+v)+z, where
%     v   .. marginal variance = (positive) smoothing width, and
%     lik .. lik function such that feval(lik{:},varargin{:}) yields a result.
% The smoothing results from a lower bound on the likelihood:
%   p(y|f) \ge exp( (b+z/ga)*f - f.^2/(2*ga) - h(ga)/2 )
function [varargout] = likVB(v, lik, varargin)
  [b,z] = feval(lik{:},varargin{1:2},[],zeros(size(v)),'infVB');
  f = varargin{3};                               % obtain location of evaluation
  sign_fmz = 2*(f-z>=0)-1;                    % strict sign mapping; sign(0) = 1
  g = sign_fmz.*sqrt((f-z).^2 + v) + z;
  varargin{3} = g;
  id = v==0 | abs(f./sqrt(v+eps))>1e10;     % correct asymptotics of f -> +/-Inf

  varargout = cell(nargout,1);              % allocate output, eval lik(..,g,..)
  [varargout{1:min(nargout,3)}] = feval(lik{:},varargin{1:3},[],'infLaplace');
  if nargout>0
    lp = varargout{1}; 
    varargout{1} = lp + b.*(f-g);
    varargout{1}(id) = lp(id);              % correct asymptotics of f -> +/-Inf
    if nargout>1                                              % first derivative
      dg_df = (abs(f-z)+eps)./(abs(g-z)+eps);    % stabilised dg/df for v=0, f=0
      dlp = varargout{2};
      varargout{2} = dlp.*dg_df + b.*(1-dg_df);
      varargout{2}(id) = dlp(id);           % correct asymptotics of f -> +/-Inf
      if nargout>2                                           % second derivative
        d2lp = varargout{3};
        g_e = g-z + sign_fmz*eps;
        v_g3  = v./(g_e.*g_e.*g_e);    % stabilised v./g.^3 to cover v=0 and f=0
        varargout{3} = (dlp-b).*v_g3 + d2lp.*dg_df.*dg_df;
        varargout{3}(id) = d2lp(id);        % correct asymptotics of f -> +/-Inf
        if nargout>3
          W = abs( (b-dlp)./(g-z+sign_fmz/1.5e8) );           % optimal sW value
          varargout{4} = sqrt(W);
          if nargout>4
            varargout{5} = b;
            if nargout>5
              varargout{6} = z;
            end
          end
        end
      end
    end
  end