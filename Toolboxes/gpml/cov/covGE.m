function [K,dK] = covGE(mode, par, hyp, varargin)

% Gamma Exponential covariance function.
% The covariance function is parameterized as:
%
% k(x,z) = exp(-r^gamma), r = maha(x,z)
%
% where maha(x,z) is a Mahalanobis distance and gamma is the shape parameter
% for the GE covariance. The hyperparameters are:
%
% hyp = [ hyp_maha
%         log(gamma/(2-gamma)) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-27.
%
% See also cov/covFunctions.m.

if nargin<1, mode = 'eye'; end, if nargin <2, par = []; end     % default values
if nargin<4, K = [covMaha(mode,par),'+1']; return, end

gamma = 2/(1+exp(-hyp(end)));
k = @(d2) exp(-d2.^(gamma/2));
dk = @(d2,k) -gamma/2*set_zero(k.*d2.^(gamma/2-1),d2==0);

if nargout==2
  [K,dKmaha,D2] = covMaha(mode,par,k,dk,hyp(1:end-1),varargin{:});
  dK = @(Q) dirder(Q,K,D2,dKmaha,gamma);
else
  K = covMaha(mode,par,k,dk,hyp(1:end-1),varargin{:});
end

function [dhyp,dx] = dirder(Q,K,D2,dKmaha,gamma)
  if nargout==1
    dhyp = dKmaha(Q);
  else
    [dhyp,dx] = dKmaha(Q);
  end
  Q = Q.*K; B = (gamma-2)*gamma/4 * D2.^(gamma/2) .* set_zero(log(D2),D2==0);
  dhyp = [dhyp; Q(:)'*B(:)];

function A = set_zero(A,I), A(I) = 0;