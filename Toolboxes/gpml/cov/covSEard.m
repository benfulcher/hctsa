function varargout = covSEard(varargin)

% Wrapper for Squared Exponential covariance function covSE.m.
%
% Squared Exponential covariance function with Automatic Relevance Detemination
% (ARD) distance measure. The covariance function is parameterized as:
%
% k(x,z) = sf^2 * exp(-(x-z)'*inv(P)*(x-z)/2)
%
% where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
% D is the dimension of the input space and sf2 is the signal variance. The
% hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          .
%         log(ell_D)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-27.
%
% See also covSE.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covSE','ard',[]},varargin{:});