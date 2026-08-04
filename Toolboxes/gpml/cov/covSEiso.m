function varargout = covSEiso(varargin)

% Wrapper for Squared Exponential covariance function covSE.m.
%
% Squared Exponential covariance function with isotropic distance measure. The
% covariance function is parameterized as:
%
% k(x,z) = sf^2 * exp(-(x-z)'*inv(P)*(x-z)/2) 
%
% where the P matrix is ell^2 times the unit matrix and sf^2 is the signal
% variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf)  ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-27.
%
% See also cov/covSE.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covSE','iso',[]},varargin{:});