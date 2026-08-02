function varargout = covRQiso(varargin)

% Wrapper for Rational Quadratic covariance function covRQ.m.
%
% Rational Quadratic covariance function with isotropic distance measure. The
% covariance function is parameterized as:
%
% k(x,z) = sf^2 * [1 + (x-z)'*inv(P)*(x-z)/(2*alpha)]^(-alpha)
%
% where the P matrix is ell^2 times the unit matrix, sf2 is the signal
% variance and alpha is the shape parameter for the RQ covariance. The
% hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf)
%         log(alpha) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-10-01.
%
% See also covRQ.m.

varargout = cell(max(1,nargout),1);
if nargin>0                                  % restore old hyper parameter order
  hyp = varargin{1};
  if numel(hyp)>2, varargin{1} = hyp([1:end-2,end,end-1]); end
end
[varargout{:}] = covScale({'covRQ','iso',[]},varargin{:});
if nargout>1                                 % restore old hyper parameter order
  o2 = varargout{2}; varargout{2} = @(Q) dirder(Q,o2);
end

function [dKdhyp,dKdx] = dirder(Q,dK)
  if nargout>1, [dKdhyp,dKdx] = dK(Q); else dKdhyp = dK(Q); end
  dKdhyp = dKdhyp([1:end-2,end,end-1]);