% Wrapper to apxSparse to remain backwards compatible.
%
% Note that covFITC is not a valid covariance function on its own right.
%
% Copyright (c) by Ed Snelson, Carl Edward Rasmussen 
%                                               and Hannes Nickisch, 2016-08-25.

function varargout = covFITC(varargin)
varargout = cell(nargout, 1); [varargout{:}] = apxSparse(varargin{:});