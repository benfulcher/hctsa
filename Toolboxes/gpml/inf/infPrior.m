function [post nlZ dnlZ] = infPrior(inf, prior, hyp, varargin)

% Inference with hyperparameter prior. Any inference method can be used.
% The prior specification is done by means of the prior parameter that has
% the same structure as hyp. Examples:
%  1) >> prior.cov{2}  = {@priorLaplace,mu,s2};
%     put a Laplace prior on the second covariance function hyperparameter
%     i.e. hyp.cov(2) ~ Laplace(mu,s2)
%  2) >> prior.mean{3} = {@priorDelta} or 
%     >> prior.mean{3} =  @priorDelta  or
%     >> prior.mean{3} = 'priorClamped'
%     fix the third mean function hyperparameter to its current value and
%     exclude it from optimisation
%     i.e. hyp.mean(3) ~ delta(hyp.mean(3))
%  3) >> prior.xu = repmat({{@priorGauss,0,1}},D*nu,1)
%     put an individual standard normal prior on every inducing input coordinate
%     for FITC inference
%     i.e. hyp.xu(i) ~ N(0,1) for all i = 1..D*nu
%  4) >> prior.multi{1} = {@priorGaussMulti,mu,s2, struct('xu',1:D*nu)}
%     put a joint Gaussian prior with mean vector mu and covariance matrix
%     s2 on the inducing input coordinates in FITC inference
%     i.e. hyp.xu ~ N(mu,s2)
%  5) >> prior.multi{1} = {@priorTMulti,mu,s2,nu, off+(1:D)}
%     put a joint Student's t prior with mean vector mu, covariance matrix s2
%     and nu degrees of freedom on the weights of the mean function e.g.
%     hyp.mean ~ T(mu,s2,nu) for meanLinear
%     Here, we do not have a struct as in 4) but the index in any2vec(hyp).
%  6) >> p = {@priorWeibull,lam,k}; g = @exp; dg = @exp; ig = @log;
%     >> prior.cov{2} = {@priorTransform,g,dg,ig,p}
%     put a univariate Weibull distribution on the exponentiated second
%     covariance function hyperparameter
%     i.e. exp(hyp.cov(2)) ~ Weibull(lam,k)
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Note that at most 1 prior is allowed per hyperparameter. We consider the
% specification of more than 1 prior for a hyperparameter as an error.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2016-10-26.
%
% See also infMethods.m, usagePrior.m, priorDistributions.m.

[post nlZ dnlZ] = inf(hyp, varargin{:});                     % perform inference
if ~isempty(prior)                % add hyperprior contributions to nlZ and dnlZ
  nam = fieldnames(prior);
  num = zeros(numel(any2vec(hyp)),1);      % number of priors per hyperparameter
  for i=1:numel(nam)     % iterate over kinds of hyperparameters cov/lik/mean/xu
    ni = nam{i};                                         % name of the ith field
    if strcmp(ni,'multi')                      % first catch multivariate priors
      p = prior.(ni);
      for j=1:numel(p)  % iterate over individual muti-hyperparameter components
        pj = p{j};
        if ~isempty(pj)          % only proceed if a nonempty prior is specified
          idx = pj{end}; pj = pj(1:end-1);            % grab index from cell end
          pjstr = pj{1}; if ~ischar(pjstr), pjstr = func2str(pjstr); end
          if numel(pjstr)<5 || ~strcmp(pjstr(end-4:end),'Multi')
            error('multivariate hyperpriors are called <Name>Multi')
          end
          if isstruct(idx)                           % massage idx into a vector
            idxs = vec2any(hyp,zeros(size(num)));             % structured index
            for nj = fieldnames(idx), idxs.(nj{1})( idx.(nj{1}) ) = 1; end
            idx = any2vec(idxs)>0;                  % linearise structured index
          else
            idxz = zeros(size(num)); idxz(idx) = 1; idx = idxz>0; % binary index
          end
          if sum(idx)<=1, error('multivariate priors need >1 hyperparams'), end
          num(idx) = num(idx)+1;                             % inc prior counter
          hypu = any2vec(hyp); dnlZ = any2vec(dnlZ);
          if     strncmp(pjstr,'priorClamped',12) || ...   % clamp deriv to zero
                 strncmp(pjstr,'priorDelta',10), dnlZ(idx) = 0;
          elseif strncmp(pjstr,'priorEqual',10)   || ...       % have same deriv
                 strncmp(pjstr,'priorSame',9),   dnlZ(idx) = sum(dnlZ(idx));
          else
            [lp,dlp] = feval(pj{:}, hypu(idx));    % evaluate prior distribution
            nlZ = nlZ-lp; dnlZ(idx) = dnlZ(idx)-dlp(:);    % update nlZ and dnlZ
          end
          dnlZ = vec2any(hyp,dnlZ);
        else
          error('multivariate priors should be non empty')
        end
      end
      continue                                       % jump to univariate priors
    end
    if ~isfield(hyp,ni), error(['unknown prior field ',ni,' in hyp']), end
    p = prior.(ni);
    if numel(p)~=numel(hyp.(ni)), error(['bad hyp/prior field length ',ni]), end
    for j=1:numel(p)         % iterate over individual hyperparameter components
      pj = p{j}; if ~iscell(pj) && ~isempty(pj), pj = {pj}; end   % enforce cell
      if ~isempty(pj)            % only proceed if a nonempty prior is specified
        num = vec2any(hyp,num); num.(ni)(j) = num.(ni)(j)+1; % inc prior counter
        num = any2vec(num);
        pj1str = pj{1}; if ~ischar(pj1str), pj1str = func2str(pj1str); end
        if strncmp(pj1str,'priorClamped',12) || strncmp(pj1str,'priorDelta',10)
          dnlZ.(ni)(j) = 0;                           % clamp derivative to zero
        else
          [lp,dlp] = feval(pj{:}, hyp.(ni)(j));    % evaluate prior distribution
          nlZ = nlZ-lp; dnlZ.(ni)(j) = dnlZ.(ni)(j)-dlp;   % update nlZ and dnlZ
        end
      end
    end
  end
  if any(num>1)                 % check for hypers with more than a single prior
    num = vec2any(hyp,num); nam = fieldnames(num);
    s = '';
    for i=1:numel(nam)
      idx = find(num.(nam{i})>1);
      for j=1:numel(idx)
        s = [s,sprintf('hyp.%s(%d) ',nam{i},idx(j))];
      end
    end
    error(['More than 1 prior specified for ',s(1:end-1),'.'])  
  end
end