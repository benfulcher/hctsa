function [post nlZ dnlZ] = infGaussLik(hyp, mean, cov, lik, x, y, opt)

% Exact inference for a GP with Gaussian likelihood.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2022-04-04.
%                                      File automatically generated using noweb.
%
% See also infMethods.m, cov/apx.m.

if nargin<7, opt = []; end                          % make sure parameter exists
if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Exact inference only possible with Gaussian likelihood');
end

[n, D] = size(x);
[m,dm] = feval(mean{:}, hyp.mean, x);           % evaluate mean vector and deriv
sn2 = exp(2*hyp.lik); W = ones(n,1)/sn2;            % noise variance of likGauss
K = apx(hyp,cov,x,opt);                        % set up covariance approximation
[ldB2,solveKiW,dW,dhyp,post.L] = K.fun(W); % obtain functionality depending on W

alpha = solveKiW(y-m);
post.alpha = K.P(alpha);                       % return the posterior parameters
post.sW = sqrt(W);                              % sqrt of noise precision vector
if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + ldB2 + n*log(2*pi*sn2)/2;    % -log marginal likelihood
  if nargout>2                                         % do we want derivatives?
    dnlZ = dhyp(alpha); dnlZ.mean = -dm(alpha);
    dnlZ.lik = -sn2*(alpha'*alpha) - 2*sum(dW)/sn2 + n;
  end
end
