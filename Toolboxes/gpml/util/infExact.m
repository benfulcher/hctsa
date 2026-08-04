% Wrapper to infGaussLik to remain backwards compatible.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-08-25.

function varargout = infExact(varargin)
varargout = cell(nargout, 1); [varargout{:}] = infGaussLik(varargin{:});