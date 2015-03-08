function K = covSEvlen(llen, hyp, x, z, i)

% Squared Exponential covariance function with spatially varying lengthscale.
% The covariance function is parameterized as:
%
% k(x,z) = sf^2 * sqrt(a/b)^D * exp(-(x-z)'*(x-z)/b) where
%          a = 2*len(x)*len(z)
%          b = len(x)^2 + len(xz)^2
%
% where len is the spatially varying lengthscale (here specified in the log
% domain), D is the dimension of the input data and sf^2 is the signal variance.
% The log-lengthscale function llen is supposed to be a valid GPML mean function
% with hyperparameters hyp_len.
%
% The hyperparameters of covSEvlen are:
%
% hyp = [ hyp_len
%         log(sf)  ]
%
% The covariance function has been introduced by Mark N. Gibbs in his 1997 PhD
% thesis and was later generalised by Paciorek&Shervish at NIPS 2004.
%
% Note that be setting len(x)=len(z)=ell for every input x and z, we
% recover covSEiso.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Hannes Nickisch, 2014-09-06.
%
% See also COVSEISO.M, COVFUNCTIONS.M.

nhl = feval(llen{:});   % number of hyperparameters reserved for the lengthscale
if nargin<3, K = [nhl,'+1']; return; end           % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x); nhl = eval(nhl);   % evaluate lengthscale hyperparameter string
hypl = hyp(1:nhl);                            % length scale function parameters
sf2 = exp(2*hyp(nhl+1));                                       % signal variance

% precompute squared distances and lengthscales
lx = exp(feval(llen{:},hypl,x));                         % evaluate lengthscales
if dg                                                               % vector kxx
  K = zeros(size(x,1),1); L2sum = 2*lx.^2; Lprod = lx.^2;
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x');    m = n;
    lz = lx;
  else                                                   % cross covariances Kxz
    K = sq_dist(x',z'); m = size(z,1);
    lz = exp(feval(llen{:},hypl,z));
  end
  L2sum = (lx.^2)*ones(1,m) + ones(n,1)*(lz.^2)';
  Lprod = lx*lz';
end

A = 2*Lprod./L2sum; B = exp(-K./L2sum);
if nargin<5                                                        % covariances
  K = sf2*A.^(D/2).*B;
else                                                               % derivatives
  if i<=nhl
    dlx = feval(llen{:},hypl,x,i).*lx;
    if dg
      dL2sum = 4*lx.*dlx; dLprod = dL2sum/2; 
    else
      if xeqz
        dlz = dlx;
      else
        dlz = feval(llen{:},hypl,z,i).*lz;
      end
      dL2sum = (2*lx.*dlx)*ones(1,m) + ones(n,1)*(2*lz.*dlz)';
      dLprod = dlx*lz'+lx*dlz';
    end
    dA = 2*dLprod./L2sum - dL2sum.*A./L2sum;
    dB = exp(-K./L2sum).*K.*(dL2sum./L2sum.^2);
    K = sf2*A.^(D/2-1).*((D/2)*dA.*B + A.*dB);
  elseif i==nhl+1
    K = 2*sf2*sqrt(2*Lprod./L2sum).^D.*exp(-K./L2sum);
  else
    error('Unknown hyperparameter')
  end
end