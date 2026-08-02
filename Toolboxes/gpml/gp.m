function [varargout] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys)
% Gaussian Process inference and prediction. The gp function provides a
% flexible framework for Bayesian inference and prediction with Gaussian
% processes for scalar targets, i.e. both regression and binary
% classification. The prior is Gaussian process, defined through specification
% of its mean and covariance function. The likelihood function is also
% specified. Both the prior and the likelihood may have hyperparameters
% associated with them.
%
% Two modes are possible: training or prediction: if no test cases are
% supplied, then the negative log marginal likelihood and its partial
% derivatives w.r.t. the hyperparameters is computed; this mode is used to fit
% the hyperparameters. If test cases are given, then the test set predictive
% probabilities are returned. Usage:
%
%   training: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, x, y);
% prediction: [ymu ys2 fmu fs2   ] = gp(hyp, inf, mean, cov, lik, x, y, xs);
%         or: [ymu ys2 fmu fs2 lp] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);
%
% where:
%
%   hyp      struct of column vectors of mean/cov/lik hyperparameters
%   inf      function specifying the inference method 
%   mean     prior mean function
%   cov      prior covariance function
%   lik      likelihood function
%   x        n by D matrix of training inputs
%   y        column vector of length n of training targets
%   xs       ns by D matrix of test inputs
%   ys       column vector of length nn of test targets
%
%   nlZ      returned value of the negative log marginal likelihood
%   dnlZ     struct of column vectors of partial derivatives of the negative
%               log marginal likelihood w.r.t. mean/cov/lik hyperparameters
%   ymu      column vector (of length ns) of predictive output means
%   ys2      column vector (of length ns) of predictive output variances
%   fmu      column vector (of length ns) of predictive latent means
%   fs2      column vector (of length ns) of predictive latent variances
%   lp       column vector (of length ns) of log predictive probabilities
%
%   post     struct representation of the (approximate) posterior
%            3rd output in training mode or 6th output in prediction mode
%            can be reused in prediction mode gp(.., cov, lik, x, post, xs,..)
%
% See also infMethods.m, meanFunctions.m, covFunctions.m, likFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2022-04-04.
%                                      File automatically generated using noweb.
if nargin<7 || nargin>9
  disp('Usage: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, x, y);')
  disp('   or: [ymu ys2 fmu fs2   ] = gp(hyp, inf, mean, cov, lik, x, y, xs);')
  disp('   or: [ymu ys2 fmu fs2 lp] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);')
  return
end

if size(y,2)>1      % deal with (independent) multivariate output y by recursing
  d = size(y,2); varargout = cell(nargout,1); out = cell(nargout,1);  % allocate
  for i=1:d
    in = {hyp, inf, mean, cov, lik, x, y(:,i)};
    if nargin>7, in = {in{:}, xs}; end
    if nargin>8, in = {in{:}, ys(:,i)}; end
    if i==1, [varargout{:}] = gp(in{:});    % perform inference for dimension ..
    else           [out{:}] = gp(in{:});             % .. number i in the output
      if nargin==7, no = 2;
        varargout{1} = varargout{1} + out{1};                          % sum nlZ
        if nargout>1                                                  % sum dnlZ
          varargout{2} = vec2any(hyp,any2vec(varargout{2})+any2vec(out{2}));
        end
      else no = 5;                              % concatenate ymu ys2 fmu fs2 lp
        for j=1:min(nargout,no), varargout{j} = [varargout{j},out{j}]; end
      end
      if nargout>no                                           % concatenate post
        if i==2, varargout{no+1} = {varargout{no+1},out{no+1}};
        else varargout{no+1} = {varargout{no+1}{:},out{no+1}}; end
      end  
    end
  end, return                                      % return to end the recursion
end

if isempty(mean), mean = {@meanZero}; end                     % set default mean
if ischar(mean) || isa(mean, 'function_handle'), mean = {mean}; end  % make cell
if isempty(cov), error('Covariance function cannot be empty'); end  % no default
if ischar(cov) || isa(cov,'function_handle'), cov  = {cov};  end     % make cell
cstr = cov{1}; if isa(cstr,'function_handle'), cstr = func2str(cstr); end
if (strcmp(cstr,'covFITC') || strcmp(cstr,'apxSparse')) && isfield(hyp,'xu')
  cov{3} = hyp.xu;                                                   %use hyp.xu
end
if isempty(inf), inf = {@infGaussLik}; end        % set default inference method
if ischar(inf),  inf = str2func(inf);  end        % convert into function handle
if ischar(inf) || isa(inf,'function_handle'), inf = {inf};  end      % make cell
istr = inf{1}; if isa(istr,'function_handle'), istr = func2str(istr); end
if strcmp(istr,'infPrior')
  istr = inf{2}; if isa(istr,'function_handle'), istr = func2str(istr); end
end
if isempty(lik), lik = {@likGauss}; end                        % set default lik
if ischar(lik) || isa(lik,'function_handle'), lik = {lik};  end      % make cell
lstr = lik{1}; if isa(lstr,'function_handle'), lstr = func2str(lstr); end

D = size(x,2);
if strncmp(cstr,'covGrid',7) || strcmp(cstr,'apxGrid') % only some inf* possible
  D = 0; xg = cov{3}; p = numel(xg); for i=1:p, D = D+size(xg{i},2); end  % dims
end

if ~isfield(hyp,'mean'), hyp.mean = []; end        % check the hyp specification
if eval(feval(mean{:})) ~= numel(hyp.mean)
  error('Number of mean function hyperparameters disagree with mean function')
end
if ~isfield(hyp,'cov'), hyp.cov = []; end
if eval(feval(cov{:})) ~= numel(hyp.cov)
  error('Number of cov function hyperparameters disagree with cov function')
end
if ~isfield(hyp,'lik'), hyp.lik = []; end
if eval(feval(lik{:})) ~= numel(hyp.lik)
  error('Number of lik function hyperparameters disagree with lik function')
end

try                                                  % call the inference method
  % issue a warning if a classification likelihood is used in conjunction with
  % labels different from +1 and -1
  if strcmp(lstr,'likErf') || strcmp(lstr,'likLogistic')
    if ~isstruct(y)
      uy = unique(y);
      if any( uy~=+1 & uy~=-1 )
        warning('You try classification with labels different from {+1,-1}')
      end
    end
  end
  if nargin>7   % compute marginal likelihood and its derivatives only if needed
    if isstruct(y)
      post = y;            % reuse a previously computed posterior approximation
    else
      post = feval(inf{:}, hyp, mean, cov, lik, x, y);
    end
  else
    if nargout<=1
      [post nlZ] = feval(inf{:}, hyp, mean, cov, lik, x, y); dnlZ = {};
    else
      [post nlZ dnlZ] = feval(inf{:}, hyp, mean, cov, lik, x, y);
    end
  end
catch
  msgstr = lasterr;
  if nargin>7, error('Inference method failed [%s]', msgstr); else 
    warning('Inference method failed [%s] .. attempting to continue',msgstr)
    varargout = {NaN, vec2any(hyp,zeros(numel(any2vec(hyp)),1))}; return % go on
  end
end

if nargin==7                                     % if no test cases are provided
  varargout = {nlZ, dnlZ, post};    % report -log marg lik, derivatives and post
else
  alpha = post.alpha; L = post.L; sW = post.sW;
  if issparse(alpha)                  % handle things for sparse representations
    nz = alpha ~= 0;                                 % determine nonzero indices
    if issparse(L), L = full(L(nz,nz)); end      % convert L and sW if necessary
    if issparse(sW), sW = full(sW(nz)); end
  else nz = true(size(alpha,1),1); end               % non-sparse representation
  if isempty(L)                       % in case L is not provided, we compute it
    K = feval(cov{:}, hyp.cov, x(nz,:));
    L = chol(eye(sum(nz))+sW*sW'.*K);
  end
  %verify whether L contains valid Cholesky decomposition or something different
  Lchol = isnumeric(L) && all(all(tril(L,-1)==0)&diag(L)'>0&isreal(diag(L))');
  ns = size(xs,1);                                       % number of data points
  if strncmp(cstr,'apxGrid',7), xs = apxGrid('idx2dat',cov{3},xs); end  % expand
  nperbatch = 1000;                       % number of data points per mini batch
  nact = 0;                       % number of already processed test data points
  ymu = zeros(ns,1); ys2 = ymu; fmu = ymu; fs2 = ymu; lp = ymu;   % allocate mem
  while nact<ns               % process minibatches of test cases to save memory
    id = (nact+1):min(nact+nperbatch,ns);               % data points to process
    kss = feval(cov{:}, hyp.cov, xs(id,:), 'diag');              % self-variance
    if strcmp(cstr,'covFITC') || strcmp(cstr,'apxSparse')    % cross-covariances
      Ks = feval(cov{:}, hyp.cov, x, xs(id,:)); Ks = Ks(nz,:); % res indep. of x
    else
      Ks = feval(cov{:}, hyp.cov, x(nz,:), xs(id,:));        % avoid computation
    end
    ms = feval(mean{:}, hyp.mean, xs(id,:));
    N = size(alpha,2);  % number of alphas (usually 1; more in case of sampling)
    Fmu = repmat(ms,1,N) + Ks'*full(alpha(nz,:));        % conditional mean fs|f
    fmu(id) = sum(Fmu,2)/N;                                   % predictive means
    if Lchol    % L contains chol decomp => use Cholesky parameters (alpha,sW,L)
      V  = L'\(repmat(sW,1,length(id)).*Ks);
      fs2(id) = kss - sum(V.*V,1)';                       % predictive variances
    else                % L is not triangular => use alternative parametrisation
      if isnumeric(L), LKs = L*Ks; else LKs = L(Ks); end    % matrix or callback
      fs2(id) = kss + sum(Ks.*LKs,1)';                    % predictive variances
    end
    fs2(id) = max(fs2(id),0);   % remove numerical noise i.e. negative variances
    Fs2 = repmat(fs2(id),1,N);     % we have multiple values in case of sampling
    if nargin<9
      [Lp, Ymu, Ys2] = feval(lik{:},hyp.lik,[],   Fmu(:),Fs2(:));
    else
      Ys = repmat(ys(id),1,N);
      [Lp, Ymu, Ys2] = feval(lik{:},hyp.lik,Ys(:),Fmu(:),Fs2(:));
    end
    lp(id)  = sum(reshape(Lp, [],N),2)/N;    % log probability; sample averaging
    ymu(id) = sum(reshape(Ymu,[],N),2)/N;          % predictive mean ys|y and ..
    ys2(id) = sum(reshape(Ys2,[],N),2)/N;                          % .. variance
    nact = id(end);          % set counter to index of last processed data point
  end
  if nargin<9
    varargout = {ymu, ys2, fmu, fs2, [], post};        % assign output arguments
  else
    varargout = {ymu, ys2, fmu, fs2, lp, post};
  end
end
