% Wrapper to apxGrid to remain backwards compatible.
%
% Note that covGrid is not a valid covariance function on its own right.
%
% Copyright (c) by Hannes Nickisch and Andrew Wilson 2016-08-25.

function varargout = covGrid(varargin)
varargout = cell(nargout, 1); [varargout{:}] = apxGrid(varargin{:});
