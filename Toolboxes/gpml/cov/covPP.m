function varargout = covPP(mode, par, v, hyp, x, varargin)

% Piecewise Polynomial covariance function with compact support, v = 0,1,2,3.
% The covariance functions are 2v times contin. diff'ble and the corresponding
% processes are hence v times  mean-square diffble. The covariance function is:
%
% k(x,z) = max(1-r,0)^(j+v) * f(r,j) with j = floor(D/2)+v+1
%
% where r is the Mahalanobis distance sqrt(maha(x,z)). The hyperparameters are:
%
% hyp = [ hyp_maha ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-05-23.
%
% See also cov/covMaha.m.

if nargin < 1, error('Mode cannot be empty.'); end                  % no default
if nargin < 2, par = []; end                                           % default
varargout = cell(max(1, nargout), 1);                  % allocate mem for output
if nargin<5, varargout{1} = covMaha(mode,par); return, end

[n,D] = size(x);
if all(v~=[0,1,2,3]), error('only 0,1,2 and 3 allowed for v'), end      % degree
j = floor(D/2)+v+1;                                                   % exponent

switch v
  case 0,  f = @(r,j) 1;
          df = @(r,j) 0;
  case 1,  f = @(r,j) 1 + (j+1)*r;
          df = @(r,j)     (j+1);
  case 2,  f = @(r,j) 1 + (j+2)*r +   (  j^2+ 4*j+ 3)/ 3*r.^2;
          df = @(r,j)     (j+2)   + 2*(  j^2+ 4*j+ 3)/ 3*r;
  case 3,  f = @(r,j) 1 + (j+3)*r +   (6*j^2+36*j+45)/15*r.^2 ...
                                + (j^3+9*j^2+23*j+15)/15*r.^3;
          df = @(r,j)     (j+3)   + 2*(6*j^2+36*j+45)/15*r    ...
                                + (j^3+9*j^2+23*j+15)/ 5*r.^2;
end
 cs = @(r,e) (r<1).*max(1-r,0).^e;
 pp = @(r,j,v,f)    cs(r,j+v  ).*  f(r,j);
dpp = @(r,j,v,f) r.*cs(r,j+v-1).*( f(r,j)*(j+v) - max(1-r,0).*df(r,j) );

k = @(d2) pp( sqrt(d2), j, v, f );
dk = @(d2,k) set_zero( -(1/2)*dpp( sqrt(d2), j, v, f )./d2 , d2==0);

[varargout{:}] = covMaha(mode, par, k, dk, hyp, x, varargin{:});
function A = set_zero(A,I), A(I) = 0;