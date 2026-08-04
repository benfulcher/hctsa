function varargout = binaryGP(hyper, approx, covfunc, lik, x, y, xstar)

% Approximate binary Gaussian Process classification. Two modes are possible:
% training or testing: if no test cases are supplied, then the approximate
% negative log marginal likelihood and its partial derivatives wrt the
% hyperparameters is computed; this mode is used to fit the hyperparameters. If
% test cases are given, then the test set predictive probabilities are
% returned. Exact inference is intractible, the function uses a specified
% approximation method (see approximations.m), flexible covariance functions
% (see covFunctions.m) and likelihood functions (see likelihoods.m).
%
% usage: [nlZ, dnlZ  ] = binaryGP(hyper, approx, covfunc, lik, x, y);
%    or: [p,mu,s2,nlZ] = binaryGP(hyper, approx, covfunc, lik, x, y, xstar);
%
% where:
%
%   hyper    is a column vector of hyperparameters
%   approx   is a function specifying an approximation method for inference 
%   covfunc  is the name of the covariance function (see below)
%   lik      is the name of the likelihood function
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of binary +1/-1 targets
%   xstar    is a nn by D matrix of test inputs
%   nlZ      is the returned value of the negative log marginal likelihood
%   dnlZ     is a (column) vector of partial derivatives of the negative
%               log marginal likelihood wrt each hyperparameter
%   p        is a (column) vector (of length nn) of predictive probabilities
%   mu       is a (column) vector (of length nn) of predictive latent means
%   s2       is a (column) vector (of length nn) of predictive latent variances
%
% The length of the vector of hyperparameters depends on the covariance
% function, as specified by the "covfunc" input to the function, specifying the
% name of a covariance function. A number of different covariance function are
% implemented, and it is not difficult to add new ones. See covFunctions.m for
% the details.
%
% The "lik" input argument specifies the name of the likelihood function (see
% likelihoods.m).
%
% The "approx" input argument to the function specifies an approximation method
% (see approximations.m). An approximation method returns a representation of
% the approximate Gaussian posterior. Usually, the approximate posterior admits
% the form N(m=K*alpha, V=inv(inv(K)+W)), where alpha is a vector and W is
% diagonal. The approximation method returns:
%
%   alpha    is a (sparse or full column vector) containing inv(K)*m, where K
%               is the prior covariance matrix and m the approx posterior mean
%   sW       is a (sparse or full column) vector containing diagonal of sqrt(W)
%               the approximate posterior covariance matrix is inv(inv(K)+W)
%   L        is a (sparse or full) matrix, L = chol(sW*K*sW+eye(n))
%
% In cases where the approximate posterior variance does not admit the form
% V=inv(inv(K)+W) with diagonal W, L contains instead -inv(K+inv(W)), and sW
% is unused.
%
% The alpha parameter may be sparse. In that case sW and L can either be sparse
% or full (retaining only the non-zero rows and columns, as indicated by the
% sparsity structure of alpha).  The L paramter is allowed to be empty, in
% which case it will be computed.
%
% The function can conveniently be used with the "minimize" function to train
% a Gaussian Process, eg:
%
% [hyper, fX, i] = minimize(hyper, 'binaryGP', length, 'approxEP', 'covSEiso', 'logistic', x, y);
%
% where "length" gives the length of the run: if it is positive, it gives the 
% maximum number of line searches, if negative its absolute gives the maximum 
% allowed number of function evaluations.
%
% Copyright (c) 2007 Carl Edward Rasmussen and Hannes Nickisch, 2010-07-21.

if nargin<6 || nargin>7
  disp('Usage: [nlZ, dnlZ  ] = binaryGP(hyper,approx,covfunc,lik,x,y);')
  disp('   or: [p,mu,s2,nlZ] = binaryGP(hyper,approx,covfunc,lik,x,y,xstar);')
  return
end

% Note, this function is just a wrapper provided for backward compatibility,
% the functionality is now provided by the more general gp function.
if strcmp(lik,'logistic')
  lik = 'likLogistic'; 
elseif strcmpi(lik,'cumgauss')
  lik = 'likErf'; 
else
  error('Allowable likelihoods: logistic and cumGauss.');
end
if strcmp(approx,'approxLA')
  inf = 'infLaplace'; 
elseif strcmp(approx,'approxEP')
  inf = 'infEP'; 
else
  error('Allowable approximations: approxLA and approxEP.');
end
mean = 'meanZero';

varargout = cell(nargout, 1);    % allocate the right number of output arguments
hyp.cov = hyper;

if nargin==6
  [varargout{:}] = gp(hyp,inf,mean,covfunc,lik,x,y);
  if nargout>1, varargout{2} = varargout{2}.cov; end
else
  [ymu,ys2,fmu,fs2] = gp(hyp,inf,mean,covfunc,lik,x,y,xstar);
  if nargout>0, varargout{1} = (1+ymu)/2; end
  if nargout>1, varargout{2} = fmu; end
  if nargout>2, varargout{3} = fs2; end
  if nargout>3, varargout{4} = gp(hyp,inf,mean,covfunc,lik,x,y); end       % nlZ
end