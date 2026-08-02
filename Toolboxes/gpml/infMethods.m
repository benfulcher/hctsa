% Inference methods: Compute the (approximate) posterior for a Gaussian process.
% Methods currently implemented include:
%
%   infGaussLik      Exact inference (only possible with Gaussian likelihood)
%   infLaplace       Laplace's Approximation
%   infEP            Expectation Propagation
%   infVB            Variational Bayes Approximation
%   infKL            Kullback-Leibler optimal Approximation
%
%   infMCMC     Markov Chain Monte Carlo and Annealed Importance Sampling
%               We offer two samplers.
%                 - hmc: Hybrid Monte Carlo
%                 - ess: Elliptical Slice Sampling
%               No derivatives w.r.t. to hyperparameters are provided.
%
%   infLOO      Leave-One-Out predictive probability and Least-Squares Approxim.
%   infPrior    Perform inference with hyperparameter prior.
%
% The interface to the approximation methods is the following:
%
%   function [post nlZ dnlZ] = inf..(hyp, mean, cov, lik, x, y)
%
% where:
%
%   hyp      is a struct of hyperparameters
%   mean     is the name of the mean function       (see meanFunctions.m)
%   cov      is the name of the covariance function (see covFunctions.m)
%   lik      is the name of the likelihood function (see likFunctions.m)
%   x        is a n by D matrix of training inputs
%   y        is a (column) vector (of size n) of targets
%
%   nlZ      is the returned value of the negative log marginal likelihood
%   dnlZ     is a (column) vector of partial derivatives of the negative
%               log marginal likelihood w.r.t. each hyperparameter
%   post     struct representation of the (approximate) posterior containing
%     alpha  is a (sparse or full column vector) containing inv(K)*(mu-m),
%               where K is the prior covariance matrix, m the prior mean,
%               and mu the approx posterior mean
%     sW     is a (sparse or full column) vector containing diagonal of sqrt(W)
%               the approximate posterior covariance matrix is inv(inv(K)+W)
%     L      is a (sparse or full) triangular matrix, L = chol(sW*K*sW+eye(n)),
%               or a full matrix, L = -inv(K+inv(W)),
%               or a function L(A) of a matrix A such that -(K+inv(W))*L(A) = A
%
% Usually, the approximate posterior to be returned admits the form
% N(mu=m+K*alpha, V=inv(inv(K)+W)), where alpha is a vector and W is diagonal.
%
% For more information on the individual approximation methods and their
% implementations, see the separate inf??.m files. See also gp.m.

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2022-04-04.
%                                      File automatically generated using noweb.

help infMethods
