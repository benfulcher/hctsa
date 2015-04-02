function K = covCos(hyp, x, z, i)

% Stationary covariance function for a sinusoid with period p in 1d:
%
% k(x,z) = sf^2*cos(2*pi*(x-z)/p)
%
% where the hyperparameters are:
%
% hyp = [ log(p)
%         log(sf) ]
%
% Note that covPeriodicNoDC converges to covCos as ell goes to infinity.
%
% Copyright (c) by James Robert Lloyd, 2013-08-05.
%
% See also COVFUNCTIONS.M, COVPERIODICNODC.M.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
if D>1, error('Covariance is defined for 1d data only.'), end
p   = exp(hyp(1));
sf2 = exp(2*hyp(2));

% precompute distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = repmat(x,1,n) - repmat(x',n,1);
  else                                                   % cross covariances Kxz
    K = repmat(x,1,size(z,1)) - repmat(z',n,1);
  end
end

K = 2*pi*K/p;
if nargin<4                                                        % covariances
    K = sf2*cos(K);
else                                                               % derivatives
  if i==1
    K = sf2*sin(K).*K;
  elseif i==2
    K = 2*sf2*cos(K);
  else
    error('Unknown hyperparameter')
  end
end