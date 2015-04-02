function K = covPPiso(v, hyp, x, z, i)

% Piecewise polynomial covariance function with compact support, v = 0,1,2,3.
% The covariance functions are 2v times contin. diff'ble and the corresponding
% processes are hence v times  mean-square diffble. The covariance function is:
%
% k(x^p,x^q) = sf^2 * max(1-r,0)^(j+v) * f(r,j) with j = floor(D/2)+v+1
%
% where r is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell^2 times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<3, K = '2'; return; end                  % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

[n,D] = size(x);
ell = exp(hyp(1));
sf2 = exp(2*hyp(2));
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
 pp = @(r,j,v,f) max(1-r,0).^(j+v).*f(r,j);
dpp = @(r,j,v,f) max(1-r,0).^(j+v-1).*r.*( (j+v)*f(r,j) - max(1-r,0).*df(r,j) );

% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sqrt( sq_dist(x'/ell) );
  else                                                   % cross covariances Kxz
    K = sqrt( sq_dist(x'/ell,z'/ell) );
  end
end

if nargin<5                                                        % covariances
  K = sf2*pp( K, j, v, f );
else                                                               % derivatives
  if i==1
    K = sf2*dpp( K, j, v, f );
  elseif i==2
    K = 2*sf2*pp( K, j, v, f );
  else
    error('Unknown hyperparameter')
  end
end