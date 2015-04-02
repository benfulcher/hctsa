% Compute infEP part of a likelihood function based on the infLaplace part using
% Gaussian-Hermite quadrature.
%
% The function is used in GLM likelihoods such as likPoisson, likGamma, likBeta
% and likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.

function varargout = lik_epquad(lik,hyp,y,mu,s2)
  n = max([length(y),length(mu),length(s2)]); on = ones(n,1);
  N = 20; [t,w] = gauher(N); oN = ones(1,N); lw = ones(n,1)*log(w');
  y = y(:).*on; mu = mu(:).*on; sig = sqrt(s2(:)).*on;            % vectors only
  [lpi,dlpi,d2lpi] = feval(lik{:},hyp,y*oN,sig*t'+mu*oN,[],'infLaplace');
  lZ = s(lpi+lw);
  dlZ = {}; d2lZ = {};
  if nargout>1                                         % 1st derivative wrt mean
    % Using p*dlp=dp, p=exp(lp), Z=sum_i wi*pi, dZ = sum_i wi*dpi we obtain
    %   dlZ = sum_i exp(lpi-lZ+lwi)*dlpi = sum_i ai*dlpi.
    a = exp(lpi - lZ*oN + lw);
    dlZ = sum(a.*dlpi,2);
    if nargout>2                                       % 2nd derivative wrt mean
      % Using d2lZ=(d2Z*Z-dZ^2)/Z^2 <=> d2Z=Z*(d2lZ+dlZ^2) and
      % d2Z = sum_i wi*d2Zi, we get d2lZ = sum_i ai*(d2lpi+dlpi^2)-dlZ^2.
      d2lZ = sum(a.*(d2lpi+dlpi.*dlpi),2) - dlZ.*dlZ;
    end
  end
  varargout = {lZ,dlZ,d2lZ};

% computes y = log( sum(exp(x),2) ), the softmax in a numerically safe way by
%  subtracting the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows
function [y,x] = s(logx)
  N = size(logx,2); max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;