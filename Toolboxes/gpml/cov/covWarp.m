function [K,dK] = covWarp(cov, p, dp, Dp, hyp, x, z)

% Apply a covariance function to p(x) rather than x i.e. warp the inputs.
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls another covariance function to do the actual work.
%
% The function computes:
%   k(x,z) = k0(p(x),p(z))
% Example:
%   k0  = {@covSEiso};
%   p = @(x) sum(2*x,2);
%   dp = @(x) 2*x;
%   Dp = 1;
%   k = {@covWarp,k0,p,dp,Dp};
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-11-14.
%
% See also covFunctions.m, covMask.m.

nh_string = feval(cov{:});    % number of hyperparameters of the full covariance
D = Dp;                                                % make variable available
if nargin<6, K = num2str(eval(nh_string)); return, end    % number of parameters
if nargin<7, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode
if numel(p) ==0,  p = @(x) x; end                          % default is identity
if numel(dp)==0, dp = @(x) 1; end                               % default is one

if eval(nh_string)~=length(hyp)                          % check hyperparameters
  error('number of hyperparameters does not match size of warped data')
end

px = p(x); if ~dg && ~xeqz, pz = p(z); else pz = z; end
if nargout>1
  [K,dK] = feval(cov{:}, hyp, px, pz);
  dK = @(Q) dirder(Q,dK,x,dp);
else
  K = feval(cov{:}, hyp, px, pz);
end

function [dhyp,dx] = dirder(Q,dK,x,dp)
  if nargout>1
    [dhyp,dx] = dK(Q); Dp = size(dx,2); dpx = dp(x);           % size(dx)=[n,Dp]
    if ndims(dpx)<=2                                          % apply chain rule
      dx = bsxfun(@times,dx,dpx);            % size(dpx)=[1,1] or [n,1] or [n,D]
    else
      dx = sum(bsxfun(@times,reshape(dx,[],1,Dp),dpx),3);   % size(dpx)=[n,D,Dp]
    end                                                         % size(dx)=[n,D]
  else
    dhyp = dK(Q);
  end