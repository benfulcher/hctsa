function varargout = covSEproj(d,varargin)

% Wrapper for Squared Exponential covariance function covSE.m.
%
% Projected squared exponential covariance. The covariance function is 
% parameterized as:
%
% k(x,z) = sf^2 * exp(-(x-z)'*inv(P)*(x-z)/2) 
%
% where the inv(P) matrix is L'*L with L of size (d,D) and sf^2 is the signal
% variance.
%
% The hyperparameters are:
% hyp = [ L_11;
%         L_21;
%          ..
%         L_d1;
%          ..
%         L_dD;
%         log(sf) ],
%
% Copyright (c) by Roman Garnett & Hannes Nickisch, 2016-04-27.
%
% See also covFunctions.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covSE','proj',d},varargin{:});