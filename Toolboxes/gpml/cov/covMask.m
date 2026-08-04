function [K,dK] = covMask(mask, cov, hyp, x, z)

% Apply a covariance function to a subset of the dimensions only. The subset can
% either be specified by a 0/1 mask by a boolean mask or by an index set.
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls another covariance function to do the actual work.
%
% The function computes:
%   k(x,z) = k0(x(mask),z(mask))
% Example:
%   k0  = {@covSEiso};
%   msk = [1,3,7];
%   k = {@covMask,msk,k0{:}};
%
% The function was suggested by Iain Murray, 2010-02-18 and is based on an
% earlier implementation of his dating back to 2009-06-16.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-11-14.
%
% See also covFunctions.m.

narg = nargin;    % make a copy to become independent of actual number of params
if iscell(mask) && numel(mask)==2  % => [K,dK] = covMask({mask, cov}, hyp, x, z)
  if narg>3, z = x; end                                  % shift parameters by 1
  if narg>2, x = hyp; end
  if narg>1, hyp = cov; end
  narg = narg+1;
  cov = mask{2}; mask = mask{1};           % split {mask, cov} into constituents
end

if ~iscell(cov), cov = {cov}; end                % properly wrap into cell array
nh_string = feval(cov{:});    % number of hyperparameters of the full covariance

if max(mask)<2 && length(mask)>1, mask = find(mask); end    % convert 1/0->index
D = length(mask);                                             % masked dimension
if narg<4, K = num2str(eval(nh_string)); return, end      % number of parameters
if narg<5, z = []; end                                     % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

if eval(nh_string)~=length(hyp)                          % check hyperparameters
  error('number of hyperparameters does not match size of masked data')
end

xm = x(:,mask); if ~dg && ~xeqz, zm = z(:,mask); else zm = z; end
if nargout>1
  [K,dK] = feval(cov{:}, hyp, xm, zm);
  dK = @(Q) dirder(Q,dK,x,mask);
else
  K = feval(cov{:}, hyp, xm, zm);
end

function [dhyp,dx] = dirder(Q,dK,x,mask)
  if nargout>1
    [dhyp,dxm] = dK(Q); n = size(x,1);
    subs = [repmat((1:n)',length(mask),1), reshape(repmat(mask(:)',n,1),[],1)];
    dx = accumarray(subs,dxm(:),size(x));
  else
    dhyp = dK(Q);
  end