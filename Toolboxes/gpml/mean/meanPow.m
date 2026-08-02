function [m,dm] = meanPow(d, mean, hyp, x)

% meanPow - compose a mean function as the power of another mean function m0.
%
% m(x) = m_0(x) ^ d
%
% If the degree d is not a strictly positive integer, we use 
% m(x) = sign(m_0(x)) * abs(m_0(x)) ^ d
% to stay within the reals
%
% The hyperparameter is:
%
% hyp = [ hyp_m0 ]
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean function to do the actual work.
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2016-05-04.
%
% See also meanFunctions.m.

if nargin<4                                        % report number of parameters
  m = feval(mean{:}); return
end

[m0,dm0] = feval(mean{:},hyp,x);                                   % evaluate m0
if d>0 && ceil(d)==d                                 % strictly positive integer
  s = 1;
else                                                       % general real number
  s = sign(m0); m0 = abs(m0); 
end
m = s.*m0.^d;                                                             % mean
dm = @(q) dm0( (d*m0.^(d-1)).*q );                      % directional derivative