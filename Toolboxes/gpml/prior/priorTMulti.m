function [lp,dlp] = priorTMulti(mu,s2,nu,x)

% Multivariate Laplace hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = ( 1 + r2/(nu-2) )^(-(nu+D)/2) / Z, where
%      Z = gamma(nu/2) * sqrt(det(pi*(nu-2)*s2)) / gamma((nu+D)/2), and
%  r2(x) = (x-mu)'*inv(s2)*(x-mu),
%
% further mu(Dx1) is the mean parameter, s2(Dx1) or s2(DxD) is the variance
% parameter, nu(1x1) > 2 is the degrees of freedom parameter and x(DxN) contains
% query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Hannes Nickisch, 2014-10-16.
%
% See also PRIORDISTRIBUTIONS.M, PRIORLAPLACE.M.

if nargin<3, error('mu, s2 and nu parameters need to be provided'), end
if ndims(mu)~=2 || size(mu,2)~=1, error('mu needs to be (Dx1)'), end
D = size(mu,1);
s2_ok = ndims(s2)==2 && all(size(s2)==[D,1] | size(s2)==[D,D]);
if ~s2_ok, error('s2 needs to be (DxD) or (Dx1)'), end
if size(s2,2)==D                                        % full multivariate case
  s = chol(s2)'; lds = sum(log(diag(s)));                 % lds = log(det(s2))/2
else                                                       % diagonal covariance
  s = sqrt(s2);  lds = sum(log(s));
end
if ~isscalar(nu), error('nu parameter needs to be a scalar'), end

if nargin<4                                             % return a random sample
  lp = randn(D,1);
  for d=1:D, lp(d) = sqrt((nu-2)/priorGamma(nu/2,2)) * lp(d); end  % unit sample
  if size(s,2)==D, lp = s*lp+mu; else lp = s.*lp+mu; end % affine transformation
  return
end
if D==1 && size(x,1)>1                            % both mu/s2 scalar => inflate
  D = size(x,1); mu = mu*ones(D,1); s = s*ones(D,1); lds = D*lds;
end
if ~(ndims(x)==2 && size(x,1)==D), error('x needs to be (Dxn)'), end

oN = ones(1,size(x,2)); oD = ones(D,1);
if size(s,2)==D
  xm = x-mu*oN; xs = s\xm;
else
  xm = x-mu*oN; xs = xm./(s*oN);
end

lZ = gammaln((nu+D)/2) - gammaln(nu/2) - D*log(pi*(nu-2))/2 - lds;
r2 = sum(xs.^2,1);
lp = -(nu+D)/2*log(1+r2/(nu-2)) + lZ;
dlp = -(nu+D)/(nu-2)*xs./(oD*(1+r2/(nu-2)));

if size(s,2)==D, dlp = s'\dlp; else dlp = dlp./(s*oN); end