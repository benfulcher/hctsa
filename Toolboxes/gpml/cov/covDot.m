function [K,dK,S] = covDot(mode, par, k, dk, hyp, x, z)

% Dot product-based covariance function. The covariance function is
% parameterized as:
%
% k(x,z) = k(s), s = dot(x,z) = x'*inv(P)*z
%
% where the matrix P is the metric. 
%
% Parameters:
% 1) mode,par:
% We offer different modes (mode) with their respective parameters (par):
% mode =   par =   inv(P) =         hyp =  
%   'eye'    []      eye(D)           []
%   'iso'    []      ell^2*eye(D)     [log(ell)]
%   'ard'    []      diag(ell.^2)     [log(ell_1); ..; log(ell_D)]
%   'proj'   d       L'*L             [L_11; L_21; ..; L_dD]
%   'fact'   d       L'*L + diag(f)   [L_11; L_21; ..; L_dD; log(f_1); ..; log(f_D)]
%
% 2) k,dk:
% The functional form of the covariance is governed by two functions:
% k:  s        -> k(x,z), s = dot(x,z) = x'*inv(P)*z
% dk: s,k(x,z) -> d k(x,z) / d s
% For example, the linear covariance uses
%   k = @(s) (c+s).^d; dk = @(s) d*(c+s).^(d-1);
% Note that not all functions k,dk give rise to a valid i.e. positive
% semidefinite covariance function k(x,z).
%
% 3) hyp,x,z:
% These input parameters follow the usual covariance function interface. For the
% composition of hyp, see 1).
%
% 4) K,dK:
% See the usual covariance function interface.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2017-09-26.
%
% See also covFunctions.m.

if nargin<1, mode = 'eye'; end, if nargin <2, par = []; end     % default values
mode_list = '''eye'', ''iso'', ''ard'', ''proj'', or ''fact''.';
switch mode
  case 'eye',  ne = '0';
  case 'iso',  ne = '1';
  case 'ard',  ne = 'D';
  case 'proj', ne = [num2str(par),'*D'];
  case 'fact', ne = [num2str(par),'*D+D'];
  otherwise,   error('Parameter mode is either %s.',mode_list)
end

if nargin<6, K = ne; return; end                   % report number of parameters
if nargin<7, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');             % sort out different modes
[n,D] = size(x); ne = eval(ne);                                     % dimensions
hyp = hyp(:);                                 % make sure hyp is a column vector
if numel(hyp)~=ne, error('Wrong number of hyperparameters'), end

switch mode                                           % mvm with metric A=inv(P)
  case 'eye',  A = @(x) x;
               dAdhyp = @(T) zeros(0,1);
  case 'iso',  iell2 = exp(-2*hyp); A = @(x) x*iell2;
               dAdhyp = @(T) -2*iell2*trace(T);
  case 'ard',  iell2 = exp(-2*hyp); A = @(x) bsxfun(@times,x,iell2'); 
               dAdhyp = @(T) -2*iell2.*diag(T);
  case 'proj', d = par; L = reshape(hyp,d,D); A = @(x) (x*L')*L;
               dAdhyp = @(T) reshape(L*(T+T'),d*D,1);
  case 'fact', d = par; L = reshape(hyp(1:d*D),d,D); f = exp(hyp(d*D+1:end));
               A = @(x) (x*L')*L + bsxfun(@times,x,f');
               dAdhyp = @(T)[reshape(L*(T+T'),d*D,1); f.*diag(T)];
end

% compute dot product
if dg                                                               % vector sxx
  Az = A(x); z = x; S = sum(x.*Az,2);
else
  if xeqz                                                 % symmetric matrix Sxx
    Az = A(x); z = x;
  else                                                         % cross terms Sxz
    Az = A(z);
  end
  S = x*Az';
end
K = k(S);                                                           % covariance
if nargout > 1
  dK = @(Q) dirder(Q,S,dk,x,Az,z,dg,xeqz,dAdhyp,mode);    % dir hyper derivative
end

function [dhyp,dx] = dirder(Q,S,dk,x,Az,z,dg,xeqz,dAdhyp,mode)
  R = dk(S).*Q; 
  switch mode
    case 'eye',  dhyp = zeros(0,1);                             % fast shortcuts
    case 'iso',  dhyp = -2*R(:)'*S(:);
    case 'ard'
      if dg
        dhyp = -2*sum(x.*bsxfun(@times,R,Az),1)';
      else
        dhyp = -2*sum(x.*(R*Az),1)';
      end
    otherwise                                              % generic computation
      if dg, T = bsxfun(@times,R,z); else T = R*z; end, T = x'*T;
      dhyp = dAdhyp(T);
  end
  if nargout > 1
    if xeqz, dx = R*Az+R'*Az; else dx = R*Az; end
  end