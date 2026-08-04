function [m,dm] = meanWSPC(d, hyp, x)

% Weighted Sum of Projected Cosines or Random Kitchen Sink features.
%
% This function represents the feature function of a zero mean GP with
% stationary covariance function. See the paper "Sparse spectrum GP regression" 
% by Lazaro-Gredilla et al., JMLR, 2010 for details.
%
% m(x) = sqrt(2/d) sum_j=1..d  a_j * cos(w_j'*x + b_j)
%
% The hyperparameter is:
% hyp = [w_1; b_1; a_1; ..; w_d; b_d; a_d]
%
% Copyright (c) by William Herlands and Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<3, m = sprintf('(D+2)*%d',d); return; end    % report number of hypers
[n,D] = size(x);
if any(length(hyp)~=eval(sprintf('(D+2)*%d',d)))
  error('Incorrect number of hyperparameters for meanRKS.')
end
hyp = reshape(hyp,D+2,d);
w = hyp(1:D,:); b = hyp(D+1,:); a = hyp(D+2,:)';    % separate hyps into w, b, a

r = bsxfun(@plus,x*w,b); cr = cos(r); m = sqrt(2/d)*cr*a;                 % mean
dm = @(q) dirder(q,x,a,r,cr,d);                         % directional derivative

function dhyp = dirder(q,x,a,r,cr,d)
  msr = -sin(r); dhyp = [x'*(msr.*(q*a')); (q'*msr).*a'; q'*cr];
  dhyp = sqrt(2/d)*dhyp(:);