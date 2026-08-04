function varargout = covPERard(varargin)

% Wrapper for PERiodic covariance function covPER.m.
%
% Periodic covariance function from an arbitrary covariance function k0 via
% embedding IR^D into IC^D.
% The covariance function is parameterized as:
%
% k(x,z) = k0(u(x),u(z)), u(x) = [sin(pi*x./p); cos(pi*x./p)]
%
% where the period p belongs to covPERiso and hyp0 belong to k0:
%
% hyp = [ log(p_1)
%         log(p_2)
%          .
%         log(p_D)
%         hyp0 ]
%
% Note that for k0 = covSEard and D = 1, a faster alternative is covPeriodic.
%
% Copyright (c) by Hannes Nickisch, 2016-04-25.
%
% See also covFunctions.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covPER('ard',varargin{:});