function K = covSEisoU(hyp, x, z, i)

% Squared Exponential covariance function with isotropic distance measure with
% unit magnitude. The covariance function is parameterized as:
%
% k(x^p,x^q) = exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2) 
%
% where the P matrix is ell^2 times the unit matrix and sf2 is the signal
% variance. The hyperparameters are:
%
% hyp = [ log(ell) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

ell = exp(hyp(1));                                 % characteristic length scale

% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x'/ell);
  else                                                   % cross covariances Kxz
    K = sq_dist(x'/ell,z'/ell);
  end
end

if nargin<4                                                        % covariances
  K = exp(-K/2);
else                                                               % derivatives
  if i==1
    K = exp(-K/2).*K;
  else
    error('Unknown hyperparameter')
  end
end