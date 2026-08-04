function varargout = covGaborard(varargin)

% Wrapper for Gabor covariance function covGabor.m.
%
% Gabor covariance function with length scales ell and periods p. The covariance
% function is parameterized as:
%
% k(x,z) = h(x-z), h(t) = exp(-sum(t.^2./(2*ell.^2)))*cos(2*pi*sum(t./p)).
%
% The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          ..
%         log(ell_D)
%         log(p_1)
%         log(p_2)
%          ..
%         log(p_D) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Note that covSM implements a weighted sum of Gabor covariance functions, but
% using an alternative (spectral) parameterization.
%
% Copyright (c) by Hannes Nickisch, 2016-04-25.
%
% See also covFunctions.m, covGabor.m, covSM.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covGabor('ard',varargin{:});