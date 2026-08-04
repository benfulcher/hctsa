function [m,dm] = meanPoly(d, hyp, x)

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
% Copyright (c) by Hannes Nickisch 2016-04-15.
%
% See also meanFunctions.m.

d = max(abs(floor(d)),1);                              % positive integer degree
if nargin<3, m = ['D*',int2str(d)]; return; end   % report number of hyperparams 

[n,D] = size(x);
a = reshape(hyp,D,d);

m = zeros(n,1);                                                % allocate memory
for j=1:d, m = m + (x.^j)*a(:,j); end                            % evaluate mean
dm = @(q) dirder(q,x,a);                                % directional derivative

function dhyp = dirder(q,x,a)
  [D,d] = size(a);
  dhyp = zeros(D*d,1); for j=1:d, dhyp((j-1)*D+(1:D)) = (x.^j)'*q(:); end