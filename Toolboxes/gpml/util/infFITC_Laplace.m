% Wrapper to infLaplace to remain backwards compatible.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-10-13.

function varargout = infFITC_Laplace(varargin)
varargout = cell(nargout, 1); [varargout{:}] = infLaplace(varargin{:});