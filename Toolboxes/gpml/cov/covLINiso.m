function varargout = covLINiso(varargin)

% Wrapper for Linear covariance function covLin.m.
%
% Linear covariance function with Automatic Relevance Determination (ARD). The
% covariance function is parameterized as:
%
% k(x,z) = x'*inv(P)*z
%
% where the P matrix is ell^2 times the unit matrix. The hyperparameters are:
%
% hyp = [ log(ell) ]
%
% Note that there is no bias term; use covConst to add a bias.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2017-09-26.
%
% See also cov/covLin.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covLIN('iso',[],varargin{:});