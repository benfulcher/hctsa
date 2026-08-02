function varargout = covSE(mode, par, varargin)
% COVSE Squared Exponential covariance function with unit amplitude.
%
% Report number of hyperparameters
%  s = COVSE (mode)
%  s = COVSE (mode, par)
%  s = COVSE (mode, par, hyp)
%
% Evaluation of k(x,x)
%  k = COVSE (mode, par, hyp, x)
%  k = COVSE (mode, par, hyp, x, [])
%
% Evaluation of diag(k(x,x))
%  k = COVSE (mode, par, hyp, x, 'diag')
%
% Evaluation of k(x,z)
%  k = COVSE (mode, par, hyp, x, z)
%
% Evaluation of function derivatives dk (w.r.t. hyp and x)
%  [k, dk] = COVSE (mode, par, hyp, x)
%  [k, dk] = COVSE (mode, par, hyp, x, [])
%  [k, dk] = COVSE (mode, par, hyp, x, 'diag')
%  [k, dk] = COVSE (mode, par, hyp, x, z)
%
% Call covFunctions to get an explanation of outputs in each mode.
%
% The covariance function is:
%
% k(x,z) = exp(-maha(x,z)/2)
%
% where maha(x,z) is a squared Mahalanobis distance. The function takes a "mode"
% parameter, which specifies precisely the Mahalanobis distance used, see
% covMaha. The function returns either the number of hyperparameters (with less
% than 3 input arguments) or it returns a covariance matrix and (optionally) a
% derivative function.
%
% See also COVMAHA, COVFUNCTIONS

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-05-23.

if nargin < 1, error('Mode cannot be empty.'); end                  % no default
if nargin < 2, par = []; end                                           % default
varargout = cell(max(1, nargout), 1);                  % allocate mem for output
if nargin < 4, varargout{1} = covMaha(mode,par); return, end

k = @(d2) exp(-d2/2); dk = @(d2,k) (-1/2)*k;         % covariance and derivative
[varargout{:}] = covMaha(mode, par, k, dk, varargin{:});
