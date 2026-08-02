function varargout = covMaterniso(d,varargin)

% Wrapper for Matern covariance function covMatern.m.
%
% Matern covariance function with nu = d/2 and isotropic distance measure. For
% d=1 the function is also known as the exponential covariance function or the 
% Ornstein-Uhlenbeck covariance in 1d. The covariance function is:
%
%   k(x,z) = sf^2 * f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t+t^2/3 for d=5.
% Here r is the distance sqrt((x-z)'*inv(P)*(x-z)), P is ell times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-27.
%
% See also covFunctions.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covMatern','iso',[],d},varargin{:});