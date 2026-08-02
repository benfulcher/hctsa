function varargout = covRQard(varargin)

% Wrapper for Rational Quadratic covariance function covRQ.m.
%
% Rational Quadratic covariance function with Automatic Relevance Determination
% (ARD) distance measure. The covariance function is parameterized as:
%
% k(x,z) = sf^2 * [1 + (x-z)'*inv(P)*(x-z)/(2*alpha)]^(-alpha)
%
% where the P matrix is diagonal with ARD parameters ell_1^2,...,ell_D^2, where
% D is the dimension of the input space, sf2 is the signal variance and alpha
% is the shape parameter for the RQ covariance. The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D)
%         log(sf)
%         log(alpha) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-10-01.
%
% See also cov/covRQ.m.

varargout = cell(max(1,nargout),1);
if nargin>0                                  % restore old hyper parameter order
  hyp = varargin{1};
  if numel(hyp)>2, varargin{1} = hyp([1:end-2,end,end-1]); end
end
[varargout{:}] = covScale({'covRQ','ard',[]},varargin{:});
if nargout>1                                 % restore old hyper parameter order
  o2 = varargout{2}; varargout{2} = @(Q) dirder(Q,o2);
end

function [dKdhyp,dKdx] = dirder(Q,dK)
  if nargout>1, [dKdhyp,dKdx] = dK(Q); else dKdhyp = dK(Q); end
  dKdhyp = dKdhyp([1:end-2,end,end-1]);