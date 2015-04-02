function K = covNNone(hyp, x, z, i)

% Neural network covariance function with a single parameter for the distance
% measure. The covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * asin(x^p'*P*x^q / sqrt[(1+x^p'*P*x^p)*(1+x^q'*P*x^q)])
%
% where the x^p and x^q vectors on the right hand side have an added extra bias
% entry with unit value. P is ell^-2 times the unit matrix and sf2 controls the
% signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sqrt(sf2) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

n = size(x,1);
ell2 = exp(2*hyp(1));
sf2 = exp(2*hyp(2));

sx = 1 + sum(x.*x,2);
if dg                                                               % vector kxx
  K = sx./(sx+ell2);
else
  if xeqz                                                 % symmetric matrix Kxx
    S = 1 + x*x';
    K = S./(sqrt(ell2+sx)*sqrt(ell2+sx)');
  else                                                   % cross covariances Kxz
    S = 1 + x*z'; sz = 1 + sum(z.*z,2);
    K = S./(sqrt(ell2+sx)*sqrt(ell2+sz)');
  end
end

if nargin<4                                                        % covariances
  K = sf2*asin(K);
else                                                               % derivatives
  if i==1                                                          % lengthscale
    if dg
      V = K;
    else
      vx = sx./(ell2+sx);
      if xeqz
        V = repmat(vx/2,1,n) + repmat(vx'/2,n,1);
      else  
        vz = sz./(ell2+sz); nz = size(z,1);
        V = repmat(vx/2,1,nz) + repmat(vz'/2,n,1);
      end
    end
    K = -2*sf2*(K-K.*V)./sqrt(1-K.*K);
  elseif i==2                                                        % magnitude
    K = 2*sf2*asin(K);
  else
    error('Unknown hyperparameter')
  end
end