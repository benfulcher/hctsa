function A = meanScale(mean, hyp, x, i)

% meanScale - compose a mean function as a scaled version of another one.
%
% m(x) = a * m_0(x)
%
% The hyperparameters are:
%
% hyp = [ a;
%         hyp_m0 ]
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean functions to do the actual work.
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2014-11-01.
%
% See also MEANFUNCTIONS.M.

if nargin<3                                        % report number of parameters
  A = [feval(mean{:}),'+1']; return
end

[n,D] = size(x);
a = hyp(1);
if nargin==3                                               % compute mean vector
  A = a*feval(mean{:},hyp(2:end),x);
else                                                 % compute derivative vector
  if i==1
    A = feval(mean{:},hyp(2:end),x);
  else
    A = a*feval(mean{:},hyp(2:end),x,i-1);
  end
end