function varargout = covLINard(varargin)

% Wrapper for Linear covariance function covLin.m.
%
% Linear covariance function with Automatic Relevance Determination (ARD). The
% covariance function is parameterized as:
%
% k(x,z) = x'*inv(P)*z
%
% where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
% D is the dimension of the input space. The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D) ]
%
% Note that there is no bias term; use covConst to add a bias.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-26.
%
% See also cov/covLin.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covLIN('ard',[],varargin{:});