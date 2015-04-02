function A = meanPoly(d, hyp, x, i)

% meanPoly - compose a mean function as a polynomial.
%
% The degree d has to be a strictly positive integer.
%
% m(x) = sum_i=1..D sum_j=1..d a_ij * x_i^j
%
% The hyperparameter is:
%
% hyp = [ a_11
%         a_21
%         ..
%         a_D1
%         a_12
%         a_22
%         ..
%         a_Dd]
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean function to do the actual work.
%
% Copyright (c) by Hannes Nickisch 2013-11-02.
%
% See also MEANFUNCTIONS.M.

d = max(abs(floor(d)),1);                              % positive integer degree
if nargin<3, A = ['D*',int2str(d)]; return; end   % report number of hyperparams 

[n,D] = size(x);
a = reshape(hyp,D,d);

A = zeros(n,1);                                                % allocate memory
if nargin==3                                               % compute mean vector
  for i=1:d, A = A + (x.^i)*a(:,i); end                          % evaluate mean
else                                                 % compute derivative vector
  if i<=d*D
    j = mod(i-1,D)+1;                                          % which dimension
    A = x(:,j).^((i-j)/D+1);                                        % derivative
  else
    A = zeros(n,1);
  end
end