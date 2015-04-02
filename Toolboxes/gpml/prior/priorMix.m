function [lp,dlp] = priorMix(w,priors,x)

% Mixture of Univariate or Multivariate hyperparameter prior distributions.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = sum_i w_i p_i(x),
%
% w(bx1) is the mixture weight vector (nonnegative, normalised to add up to 1),
% priors(bx1) is a cell array containing the priors to be mixed and x(1xN)
% contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also PRIORDISTRIBUTIONS.M.

if nargin<2, error('specification of w and priors required'), end
n = numel(w); w = w(:);
if abs(sum(w)-1)>1e-9, error('w needs to sum up to 1'), end
if any(w<0), error('w needs to be nonnegative'), end
for i=1:n, if ~iscell(priors{i}), priors{i} = {priors{i}}; end, end     % sample
if nargin<3, i = 1+nnz(rand>cumsum(w)); lp = feval(priors{i}{:}); return, end

lpi = zeros(numel(x),n); dlpi = zeros(numel(x),n);
for i=1:n, [lpi(:,i),dlpi(:,i)] = feval(priors{i}{:},x(:)); end

mx = max(lpi,[],2); lp = log(exp(lpi-mx*ones(1,n))*w) + mx;
pi = exp(lpi);     dlp = ((pi.*dlpi)*w)./(pi*w);