function K = covConst(hyp, x, z, i)

% covariance function for a constant function. The covariance function is
% parameterized as:
%
% k(x^p,x^q) = s2;
%
% The scalar hyperparameter is:
%
% hyp = [ log(sqrt(s2)) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '1'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

s2 = exp(2*hyp);                                                            % s2
n = size(x,1);

if dg                                                               % vector kxx
  K = s2*ones(n,1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = s2*ones(n);
  else                                                   % cross covariances Kxz
    K = s2*ones(n,size(z,1));
  end
end

if nargin>3                                                        % derivatives
  if i==1
    K = 2*K;
  else
    error('Unknown hyperparameter')
  end
end