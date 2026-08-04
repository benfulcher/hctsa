function varargout = covSEvlen(llen,varargin)

% Wrapper for Squared Exponential covariance function covSE.m.
%
% Squared Exponential covariance function with spatially varying lengthscale.
% The covariance function is parameterized as:
%
% k(x,z) = sf^2 * sqrt(a/b)^D * exp(-(x-z)'*(x-z)/b) where
%          a = 2*len(x)*len(z)
%          b = len(x)^2 + len(z)^2
%
% where len is the spatially varying lengthscale (here specified in the log
% domain), D is the dimension of the input data and sf^2 is the signal variance.
% The log-lengthscale function llen is supposed to be a valid GPML mean function
% with hyperparameters hyp_len.
%
% The hyperparameters of covSEvlen are:
%
% hyp = [ hyp_len
%         log(sf)  ]
%
% The covariance function has been introduced by Mark N. Gibbs in his 1997 PhD
% thesis and was later generalised by Paciorek&Schervish at NIPS 2004.
%
% Note that by setting len(x)=len(z)=ell for every input x and z, we
% recover covSEiso.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2016-05-04.
%
% See also cov/covSEiso.m, covFunctions.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covSE','vlen',llen},varargin{:});