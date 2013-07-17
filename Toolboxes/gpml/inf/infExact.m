function [post nlZ dnlZ] = infExact(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood. Compute a parametrization
% of the posterior, the negative log marginal likelihood and its derivatives
% w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-01-21
%
% See also INFMETHODS.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end
 
[n, D] = size(x);
K = feval(cov{:}, hyp.cov, x);                      % evaluate covariance matrix
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
L = chol(K/sn2+eye(n));               % Cholesky factor of covariance with noise
alpha = solve_chol(L,y-m)/sn2;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L  = L;                                        % L = chol(eye(n)+sW*sW'.*K)

if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + sum(log(diag(L))) + n*log(2*pi*sn2)/2;  % -log marg lik
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    Q = solve_chol(L,eye(n))/sn2 - alpha*alpha';    % precompute for convenience
    for i = 1:numel(hyp.cov)
      dnlZ.cov(i) = sum(sum(Q.*feval(cov{:}, hyp.cov, x, [], i)))/2;
    end
    dnlZ.lik = sn2*trace(Q);
    for i = 1:numel(hyp.mean), 
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*alpha;
    end
  end
end
