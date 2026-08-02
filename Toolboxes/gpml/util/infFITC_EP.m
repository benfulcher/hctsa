% Wrapper to infEP to remain backwards compatible.
%
% FITC-EP approximation to the posterior Gaussian process. The function is
% equivalent to infEP with the covariance function:
%   Kt = Q + G; G = diag(g); g = diag(K-Q);  Q = Ku'*inv(Kuu + snu2*eye(nu))*Ku;
% where Ku and Kuu are covariances w.r.t. to inducing inputs xu and
% snu2 = mean(diag(Kuu))/1e6 is the noise of the inducing inputs which we
% fixe to 0.1% of signal standard deviation.
% For details, see The Generalized FITC Approximation, Andrew Naish-Guzman and
%                  Sean Holden, NIPS, 2007.
%
% See also infMethods.m, APXSPARSE.M, APX.M.
%
% Copyright (c) by Hannes Nickisch, 2016-10-13.

function varargout = infFITC_EP(varargin)
varargout = cell(nargout, 1); [varargout{:}] = infEP(varargin{:});
