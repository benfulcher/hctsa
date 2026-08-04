function [K,dK] = covPER(mode, cov, hyp, x, z)

% Periodic covariance function from an arbitrary covariance function k0 via
% embedding IR^D into IC^D.
% The covariance function is parameterized as:
%
% k(x,z) = k0(u(x),u(z)), u(x) = [sin(pi*x/p); cos(pi*x/p)]
%
% where the period p belongs to covPER and hyp0 belongs to k0:
%
% hyp = [ hyp_p
%         hyp0 ]
%
% We offer three different modes:
%   'eye':  p =            ones(D,1); hyp_p = [];
%   'iso':  p = exp(hyp_p)*ones(D,1); hyp_p = [log(p)];
%   'ard':  p = exp(hyp_p);           hyp_p = [log(p_1); ..; log(p_D)];
%
% Note that for k0 = covSEiso and D = 1, a faster alternative is covPeriodic.
%
% Copyright (c) by Hannes Nickisch, 2016-04-25.
%
% See also covFunctions.m.

if nargin<2, error('We require a mode and a base covariance k0.'), end
if     isequal(mode,'ard'), ne = 'D';
elseif isequal(mode,'iso'), ne = '1';
elseif isequal(mode,'eye'), ne = '0';
else error('Parameter mode is either ''eye'', ''iso'' or ''ard''.'), end

if nargin<4                                        % report number of parameters
    % occurences of D replaced by 2*D as the size of u(x) is twice the size of x
  K = ['(',ne,'+',strrep(feval(cov{:}),'D','2*D'),')']; return
end
if nargin<5, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x); ne = eval(ne);
p = exp(hyp(1:ne)); if numel(p)==0, p = 1; end, p = p.*ones(D,1);   % period par
ux = u(x,p); uz = u(z,p);                   % apply the embedding u:IR^D->IR^2*D

[K,dK0] = feval(cov{:},hyp(ne+1:end),ux,uz);               % covariance function
if nargout > 1                                          % directional derivative
  dK = @(Q) dirder(Q,dK0,p,mode,dg,xeqz,x,ux,z,uz,cov,hyp(ne+1:end));
end

function [dhyp,dx] = dirder(Q,dK0,p,mode,dg,xeqz,x,ux,z,uz,cov,hypc)
  if dg                               % dx not required for dg so we assume zero
    dhyp = dK0(Q); dux = zeros(size(ux));          % only correct for stationary
  else
    [dhyp,dux] = dK0(Q);
  end
  if isequal(mode,'eye')
    dp = zeros(0,1); if nargout > 1, dx = duxdx(dux,ux,p); end
  else
    dx = duxdx(dux,ux,p); dp = -sum(dx.*x,1);
    if ~xeqz && ~dg
      [junk,dK0t] = feval(cov{:},hypc,uz,ux); [junk,duz] = dK0t(Q');
      dp = dp - sum(duxdx(duz,uz,p).*z,1);
    end
    if isequal(mode,'iso'), dp = sum(dp); end
  end
  dhyp = [dp(:); dhyp];

function ux = u(x,p)                        % apply the embedding u:IR^D->IR^2*D
  if numel(x)==0 || ischar(x)
    ux = x;
  else
    ux = 2*pi*bsxfun(@times,x,1./p(:)'); ux = [sin(ux), cos(ux)];
  end

function dx = duxdx(du,ux,p)
  D = size(ux,2)/2;
  dx =  bsxfun(@times, du(:,1:D).*ux(:,D+1:2*D), 2*pi./p(:)') ...
       -bsxfun(@times, du(:,D+1:2*D).*ux(:,1:D), 2*pi./p(:)');