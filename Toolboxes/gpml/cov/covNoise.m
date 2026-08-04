function varargout = covNoise(varargin)

% Wrapper for unit "white noise" covariance function covEye.m.
%
% Independent covariance function, i.e. "white noise".
% The covariance function is specified as:
%
% k(x^p,x^q) = sf^2 * \delta(p,q)
%
% \delta(p,q) is a Kronecker delta function which is 1 iff p=q and zero
% otherwise in mode 1).
% In cross covariance mode 2) two data points x_p and z_q are considered equal
% if their difference norm |x_p-z_q| is less than eps, the machine precision.
% The hyperparameters are:
%
% hyp = [ log(sf) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-19.
%
% See also cov/covEye.m.

varargout = cell(max(nargout,1),1);
[varargout{:}] = covScale({'covEye'},varargin{:});