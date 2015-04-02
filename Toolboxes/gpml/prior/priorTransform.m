function [lp,dlp] = priorTransform(g,dg,ig,prior,x)

% Transformed univariate hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%  prior(y), where y = g(x),
%
% g is the transformation function (monotonically increasing and invertable),
% dg is its derivative, ig is its inverse, prior is the distribution to be
% transformed and x(1xN) contains query hyperparameters for prior evaluation.
% Note that p(x) = prior(g(x)) is not necessarily normalised w.r.t. x.
%
% For example, to put a nonnegative gamma prior on the lenthscale
% parameter ell of covSEiso, one can use 
%   {@priorTransform,@exp,@exp,@log,{@priorGamma,k,t}}
% since ell is represented in the log domain.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-08.
%
% See also PRIORDISTRIBUTIONS.M.

if nargin<4, error('g, dg, ig and prior parameters need to be provided'), end
if nargin<5, lp = ig(feval(prior{:})); return, end      % apply inverse sampling

if ~iscell(prior), prior = {prior}; end
[lp,dlp] = feval(prior{:},g(x)); dlp = dlp.*dg(x);