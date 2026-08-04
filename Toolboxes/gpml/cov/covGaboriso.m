function varargout = covGaboriso(varargin)

% Wrapper for Gabor covariance function covGabor.m.
%
% Gabor covariance function with length scale ell and period p. The 
% covariance function is parameterized as:
%
% k(x,z) = h(x-z) with h(t) = exp(-t'*t/(2*ell^2))*cos(2*pi*sum(t)/p).
%
% The hyperparameters are:
%
% hyp = [ log(ell)
%         log(p)   ]
%
% Note that covSM implements a weighted sum of Gabor covariance functions, but
% using an alternative (spectral) parameterization.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2016-04-25.
%
% See also covFunctions.m, cov/covGabor.m, covSM.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covGabor('iso',varargin{:});