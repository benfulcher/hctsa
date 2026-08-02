function [K,dK] = covRQ(mode, par, hyp, varargin)

% Rational Quadratic covariance function.
% The covariance function is parameterized as:
%
% k(x,z) = [1 + maha(x,z)/(2*alpha)]^(-alpha)
%
% where maha(x,z) is a Mahalanobis distance and alpha is the shape parameter
% for the RQ covariance. The hyperparameters are:
%
% hyp = [ hyp_maha
%         log(alpha) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-05-23.
%
% See also covFunctions.m.

if nargin < 1, error('Mode cannot be empty.'); end                  % no default
if nargin < 2, par = []; end                                           % default
if nargin<4, K = [covMaha(mode,par),'+1']; return, end

alpha = exp(hyp(end));
k = @(d2) (1+0.5*d2/alpha).^(-alpha); dk = @(d2,k) -k./(2+d2/alpha);

if nargout==2
  [K,dKmaha,D2] = covMaha(mode,par,k,dk,hyp(1:end-1),varargin{:});
  dK = @(Q) dirder(Q,K,D2,dKmaha,alpha);
else
  K = covMaha(mode,par,k,dk,hyp(1:end-1),varargin{:});
end

function [dhyp,dx] = dirder(Q,K,D2,dKmaha,alpha)
  if nargout==1
    dhyp = dKmaha(Q);
  else
    [dhyp,dx] = dKmaha(Q);
  end
  Q = Q.*K; B = 1+0.5*D2/alpha;
  dhyp = [dhyp; sum(sum( Q.*(0.5*D2./B-alpha*log(B)) ))];