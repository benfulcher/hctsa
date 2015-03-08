function K = covSEfact(d, hyp, x, z, i)

% Factor analysis squared exponential covariance.
%
% This provides a GPML-compatible covariance function implementing the
% The covariance function is parameterized as:
%
%   k(x,z) = sf^2 * cosSEiso(x*L',z*L'),
%
% where x and z are the inputs and, L is a (dxD) embedding matrix, sf is the
% signal standard deviation.
%
% The hyperparameters are:
% hyp = [ hypL
%         log(sf) ],
%
% where sf is the signal standard deviation and hypL is the vectorized L matrix.
% The conversion between hypL and L is achieved by the following code:
%
% conversion L -> hypL
% if d==1, diagL = L(1); else diagL = diag(L); end
% L(1:d+1:d*d) = log(diagL); hyp = L(triu(true(d,D)));
%
% conversion hypL -> L
% L = zeros(d,D); L(triu(true(d,D))) = hyp(:);
% if d==1, diagL = L(1); else diagL = diag(L); end
% L(1:d+1:d*d) = exp(diagL);
%
% hypL = [ log(L_11)
%              L_21
%          log(L_22)
%              L_31
%              L_32
%              ..
%          log(L_dd)
%              ..
%              L_dD]
%
% Copyright (c) by Roman Garnett & Hannes Nickisch, 2014-08-15.
%
% See also COVFUNCTIONS.M.

if nargin==0, error('d must be specified.'), end
nh = sprintf('(D*%d - %d*(%d-1)/2 + 1)',d,d,d);    % number of hyperparam string
if nargin<3, K = nh; return; end              % report number of hyperparameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode
if xeqz, z = x; end                                % make sure we have a valid z
[n,D] = size(x);                                                % dimensionality
if d>D, error('We need d<=D.'), end
nh = eval(nh);
sf = exp(hyp(nh)); hyp = hyp(1:nh-1);            % bring hypers in correct shape

L = zeros(d,D); L(triu(true(d,D))) = hyp(:);                  % embedding matrix
if d==1, diagL = L(1); else diagL = diag(L); end    % properly handle limit case
L(1:d+1:d*d) = exp(diagL);

xl = x*L';                                                % transform input data
if dg, zl = 'diag';
else
  if xeqz, zl = xl; else zl = z*L'; end
end
if nargin<5                   % covariance, call covSEiso on the embedded points
  K = covSEiso([0;log(sf)],xl,zl);
else
  if i==nh                     % derivative w.r.t. log signal standard deviation
    K = covSEiso([0;log(sf)],xl,zl,2);
  else
    td = d*(d+1)/2;                                     % determine indices in L
    if i>td
      col = ceil((i-td)/d);
      col = col+ceil((sqrt(8*(i-d*col)+1)-1)/2);
      row = mod(i-d*(d+1)/2-1,d)+1;
    else
      col = ceil((sqrt(8*i+1)-1)/2); row = i-col*(col-1)/2;
    end
    if dg
      K = zeros(n,1);
    else
      K = covSEiso([0;log(sf)],xl,zl);       % derivative w.r.t. log lengthscale
      m = size(z,1);
      F = x(:,col)*ones(1,m)-ones(n,1)*z(:,col)';
      K = -K.*F.*(xl(:,row)*ones(1,m)-ones(n,1)*zl(:,row)');
      if row==col, K = L(row,col)*K; end         % diagonal is in the log domain
    end
  end
end
