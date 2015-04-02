function K = covNoise(hyp, x, z, i)

% Independent covariance function, i.e. "white noise", with specified variance.
% The covariance function is specified as:
%
% k(x^p,x^q) = s2 * \delta(p,q)
%
% where s2 is the noise variance and \delta(p,q) is a Kronecker delta function
% which is 1 iff p=q and zero otherwise in mode 1).
% In cross covariance mode 2) two data points x_p and z_q are considered equal
% if their difference norm |x_p-z_q| is less than eps, the machine precision.
% The hyperparameter is:
%
% hyp = [ log(sqrt(s2)) ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-05-15.
%
% See also COVFUNCTIONS.M.

tol = eps;   % threshold on the norm when two vectors are considered to be equal
if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
dg = strcmp(z,'diag');                                          % determine mode

n = size(x,1);
s2 = exp(2*hyp);                                                % noise variance

% precompute raw
if dg                                                               % vector kxx
  K = ones(n,1);
else
  if numel(z)==0                                          % symmetric matrix Kxx
    K = eye(n);
  else                                                   % cross covariances Kxz
    K = double(sq_dist(x',z')<tol*tol);
  end
end

if nargin<4                                                        % covariances
  K = s2*K;
else                                                               % derivatives
  if i==1
    K = 2*s2*K;
  else
    error('Unknown hyperparameter')
  end
end
