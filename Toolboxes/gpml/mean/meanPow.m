function A = meanPow(mean, hyp, x, i)

% meanPow - compose a mean function as the power of another mean function.
%
% The degree d has to be a strictly positive integer.
%
% m(x) = m_0(x) ^ d
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean function to do the actual work.
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-06-18.
%
% See also MEANFUNCTIONS.M.

d = mean{1}; d = abs(floor(d)); d = max(d,1);          % positive integer degree
mean = mean{2}; if ~iscell(mean), mean = {mean}; end
if nargin<3                                        % report number of parameters
  A = feval(mean{:}); return
end

[n,D] = size(x);
if nargin==3                                               % compute mean vector
  A = feval(mean{:},hyp,x).^d;
else                                                 % compute derivative vector
  A = ( d*feval(mean{:},hyp,x).^(d-1) ).* feval(mean{:},hyp,x,i);
end