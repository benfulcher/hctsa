function [m,dm] = meanScale(mean, hyp, x)

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
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2016-04-15.
%
% See also meanFunctions.m.

if nargin<3                                        % report number of parameters
  m = [feval(mean{:}),'+1']; return
end

a = hyp(1);
[m0,dm0] = feval(mean{:},hyp(2:end),x);                            % evaluate m0
m = a*m0;                                                                 % mean
dm = @(q) [m0'*q(:); a*dm0(q)];                         % directional derivative
