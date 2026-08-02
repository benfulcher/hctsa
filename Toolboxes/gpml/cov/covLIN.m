function [K,dK] = covLIN(mode,par,hyp,x,z)

% Linear covariance function.
% The covariance function is parameterized as:
%
% k(x,z) = dot(x,z)
%
% where dot(x,z) is a dot product. The hyperparameters are:
%
% hyp = [ hyp_dot ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2017-09-26.
%
% See also cov/covDot.m.

if nargin<2, mode = 'eye'; par = []; end, narg = nargin;          % default mode
if ~ischar(mode)                      % make compatible to old interface version
  if nargin>2, z = hyp; end
  if nargin>1, x = par; end
  if nargin>0, hyp = mode; end
  mode = 'eye'; narg = narg+2;
end
if narg<4, K = covDot(mode,par); return, end
if narg<5, z = []; end                                     % make sure, z exists

k = @(s) s; dk = @(s) ones(size(s));

if nargout > 1
  [K,dK] = covDot(mode,par,k,dk,hyp,x,z);
else
  K = covDot(mode,par,k,dk,hyp,x,z);
end