function varargout = covConst(varargin)

% Wrapper for unit constant covariance function covOne.m.
%
% Covariance function for a constant function. The covariance function is
% parameterized as:
%
% k(x,z) = sf^2
%
% The scalar hyperparameter is:
%
% hyp = [ log(sf) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-19.
%
% See also covFunctions.m, covOne.m.

varargout = cell(max(nargout,1),1);
[varargout{:}] = covScale({'covOne'},varargin{:});