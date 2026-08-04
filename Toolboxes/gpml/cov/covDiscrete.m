function [K,dK] = covDiscrete(s, hyp, x, z)

% Covariance function for discrete inputs. Given a function defined on the
% integers 1,2,3,..,s, the covariance function is parameterized as:
%
% k(x,z) = K_{xz},
%
% where K is a matrix of size (s x s).
%
% This implementation assumes that the inputs x and z are given as integers
% between 1 and s, which simply index the matrix K.
%
% The hyperparameters specify the upper-triangular part of the Cholesky factor
% of K, where the (positive) diagonal elements are specified by their logarithm.
% If L = chol(K), so K = L'*L, then the hyperparameters are:
%
% hyp = [ log(L_11)
%             L_21
%         log(L_22)
%             L_31
%             L_32
%         ..
%         log(L_ss) ]
%
% The hyperparameters hyp can be generated from K using:
% L = chol(K); L(1:(s+1):end) = log(diag(L));
% hyp = L(triu(true(s)));
%
% The covariance matrix K is obtained from the hyperparameters hyp by:
% L = zeros(s); L(triu(true(s))) = hyp(:);
% L(1:(s+1):end) = exp(diag(L)); K = L'*L;
%
% This parametrization allows unconstrained optimization of K.
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Roman Garnett, 2016-04-18.
%
% See also meanDiscrete.m, covFunctions.m.

if nargin==0, error('s must be specified.'), end           % check for dimension
if nargin<3, K = num2str(s*(s+1)/2); return; end   % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode
if xeqz, z = x; end                                % make sure we have a valid z

L = zeros(s); L(triu(true(s))) = hyp(:);                 % build Cholesky factor
L(1:(s+1):end) = exp(diag(L));

A = L'*L;   % A is a placeholder for K to avoid a name clash with the return arg
if dg
  K = A(sub2ind(size(A),fix(x),fix(x)));
else
  K = A(fix(x),fix(z));
end

if nargout > 1
  dK = @(Q) dirder(Q,s,L,x,z,dg);                 % directional hyper derivative
end

function [dhyp,dx] = dirder(Q,s,L,x,z,dg)
  nx = numel(x); Ix = sparse(fix(x),1:nx,ones(nx,1),s,nx);      % K==Ix'*L'*L*Iz
  if dg
    B = Ix*diag(Q)*Ix';
  else
    nz = numel(z); Iz = sparse(fix(z),1:nz,ones(nz,1),s,nz);
    B = Iz*Q'*Ix';
  end
  B = L*(B+B'); B(1:(s+1):end) = diag(B).*diag(L);         % fix exp on diagonal
  dhyp = B(triu(true(s)));
  if nargout > 1, dx = zeros(size(x)); end