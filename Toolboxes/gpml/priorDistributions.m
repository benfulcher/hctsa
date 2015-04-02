% prior distributions to be used for hyperparameters of Gaussian processes
% using infPrior.
% There are two different kinds of prior distributions: simple and composite:
%
% simple prior distributions:
%
%   priorGauss           - univariate Gaussian
%   priorLaplace         - univariate Laplace
%   priorT               - univariate Student's t
%
%   priorSmoothBox1      - univariate interval (linear decay in log domain)
%   priorSmoothBox2      - univariate interval (quadr. decay in log domain)
%
%   priorGamma           - univariate Gamma, IR+
%   priorWeibull         - univariate Weibull, IR+
%   priorInvGauss        - univariate Inverse Gaussian, IR+
%   priorLogNormal       - univariate Log-normal, IR+
%
%   priorClamped or      - fix hyperparameter to its current value by setting
%   priorDelta          derivatives to zero, no effect on marginal likelihood
%
%   priorGaussMulti      - multivariate Gauss
%   priorLaplaceMulti    - multivariate Laplace
%   priorTMulti          - multivariate Student's t
%
%   priorClampedMulti or - fix hyperparameter to its current value by setting
%   priorDeltaMulti     derivatives to zero, no effect on marginal likelihood
%
% composite prior distributions (see explanation at the bottom):
%
%   priorMix             - nonnegative mixture of priors
%   priorTransform       - prior on g(t) rather than t
%
% Naming convention: all prior distributions are named "prior/prior*.m".
%
%
% 1) With only a fixed input arguments:
%
%    r = priorNAME(par1,par2,parN)
%
% The function returns a random sample from the distribution for e.g.
% random restarts, simulations or optimisation initialisation.
%
% 2) With one additional input arguments:
%
%    [lp,dlp] = priorNAME(par1,par2,parN, t)
%
% The function returns the log density at location t along with its first
% derivative.
%
% See also doc/usagePrior.m, inf/infPrior.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
