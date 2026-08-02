function [K,dK] = covGabor(mode, hyp, x, z)

% Gabor covariance function with length scale ell and period p. The 
% covariance function is parameterized as:
%
% k(x,z) = h(x-z), h(t) = exp(-sum(t.^2./(2*ell.^2)))*cos(2*pi*sum(t./p)).
%
% The hyperparameters are:
%
% hyp = [ hyp_ell
%         hyp_p    ]
%
% We offer three different modes:
%   'eye':  ell =              ones(D,1);   p =            ones(D,1);
%   'iso':  ell = exp(hyp_ell)*ones(D,1);   p = exp(hyp_p)*ones(D,1);
%   'ard':  ell = exp(hyp_ell)          ;   p = exp(hyp_p)          ;
%
% Note that covSM implements a weighted sum of Gabor covariance functions, but
% using an alternative (spectral) parameterization.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2016-05-03.
%
% See also covFunctions.m, cov/covGaborard.m, cov/covSM.m.

if nargin<1, error('We require a mode.'), end
if     isequal(mode,'ard'), np = 'D';
elseif isequal(mode,'iso'), np = '1';
elseif isequal(mode,'eye'), np = '0';
else error('Parameter mode is either ''eye'', ''iso'' or ''ard''.'), end

if nargin<3, K = ['(2*',np,')']; return, end       % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x); np = eval(np);                                 % dimensionality
p = exp(hyp(np+1:2*np)); if numel(p)==0, p = 1; end, p = p.*ones(D,1);  % period
[Kse,dKse] = covSE(mode,[],hyp(1:np),x,z); % squared exponential base covariance

Dp = zeros(size(Kse));                               % init sum(t)/p computation
if ~dg
  if xeqz                                                 % symmetric matrix Kxx
    for d=1:D, Dp = Dp + bsxfun(@minus,x(:,d),x(:,d)')/p(d); end
  else                                                   % cross covariances Kxz
    for d=1:D, Dp = Dp + bsxfun(@minus,x(:,d),z(:,d)')/p(d); end
  end
end
C = cos(2*pi*Dp); K = Kse .* C;                                     % covariance
if nargout > 1                                          % directional derivative
  dK = @(Q) dirder(Q,C,Dp,Kse,dKse,np,dg,xeqz,x,z,D,p,mode);
end

function [dhyp,dx] = dirder(Q,C,Dp,Kse,dKse,np,dg,xeqz,x,z,D,p,mode)
  dhyp = zeros(2*np,1);                                        % allocate memory
  if nargout > 1
    [dhyp(1:np),dx] = dKse(Q.*C);
  else
    dhyp(1:np) = dKse(Q.*C);
  end
  if ~dg
    Q = Q .* Kse .* -sin(2*pi*Dp)*2*pi;
    if isequal(mode,'ard')
      for d=1:D
        if xeqz
          Dd = bsxfun(@minus,x(:,d),x(:,d)');
        else
          Dd = bsxfun(@minus,x(:,d),z(:,d)');
        end
        dhyp(np+d) = -(Q(:)'*Dd(:))/p(d);
      end
    elseif isequal(mode,'iso')
      dhyp(np+1) = -(Q(:)'*Dp(:));
    end
    if nargout > 1
      if xeqz
        dx = dx + (sum(Q,2)-sum(Q,1)')*(1./p)';
      else
        dx = dx + sum(Q,2)*(1./p)';
      end
    end
  end