function varargout = covMaternard(d,varargin)

% Wrapper for Matern covariance function covMatern.m.
%
% Matern covariance function with nu = d/2 and with Automatic Relevance
% Determination (ARD) distance measure. For d=1 the function is also known as
% the exponential covariance function or the Ornstein-Uhlenbeck covariance 
% in 1d. The covariance function is:
%
%   k(x,z) = f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t+t^2/3 for d=5.
% Here r is the distance sqrt((x-z)'*inv(P)*(x-z)), where the P matrix
% is diagonal with ARD parameters ell_1^2,...,ell_D^2, where D is the dimension
% of the input space and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D)
%         log(sf) ]
%
% Copyright (c) by Hannes Nickisch, 2016-04-17.
%
% See also covFunctions.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covMatern','ard',[],d},varargin{:});