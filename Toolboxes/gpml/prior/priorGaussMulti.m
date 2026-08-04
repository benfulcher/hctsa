function [lp,dlp] = priorGaussMulti(mu,s2,x)

% Multivariate Gaussian hyperparameter prior distribution.
% Compute log-likelihood and its derivative or draw a random sample.
% The prior distribution is parameterized as:
%
%   p(x) = exp(-r2/2) / sqrt(det(2*pi*s2)) where r2(x) = (x-mu)'*inv(s2)*(x-mu),
%
% mu(Dx1) is the mean parameter, s2(Dx1) or s2(DxD) is the variance parameter
% and x(DxN) contains query hyperparameters for prior evaluation.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-09-12.
%
% See also priorDistributions.m, prior/priorGauss.m.

if nargin<2, error('mu and s2 parameters need to be provided'), end
if ndims(mu)~=2 || size(mu,2)~=1, error('mu needs to be (Dx1)'), end
D = size(mu,1);
s2_ok = ndims(s2)==2 && all(size(s2)==[D,1] | size(s2)==[D,D]);
if ~s2_ok, error('s2 needs to be (DxD) or (Dx1)'), end
if size(s2,2)==D                                        % full multivariate case
  s = chol(s2)'; lds = sum(log(diag(s)));                 % lds = log(det(s2))/2
else                                                       % diagonal covariance
  s = sqrt(s2);  lds = sum(log(s));
end
if nargin<3                                             % return a random sample
  lp = randn(D,1);                                                 % unit sample
  if size(s,2)==D, lp = s*lp+mu; else lp = s.*lp+mu; end % affine transformation
  return
end
if D==1 && size(x,1)>1                            % both mu/s2 scalar => inflate
  D = size(x,1); mu = mu*ones(D,1); s = s*ones(D,1); lds = D*lds;
end
if ~(ndims(x)==2 && size(x,1)==D), error('x needs to be (Dxn)'), end

oN = ones(1,size(x,2));
if size(s,2)==D
  xs = s\(x-mu*oN);       dlp = -s'\xs;
else
  xs = (x-mu*oN)./(s*oN); dlp = -xs./(s*oN);
end
lp = -sum(xs.^2,1)/2 - D*log(2*pi)/2 - lds;