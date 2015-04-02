% Compute y = log( sum(exp(x),2) ), the softmax in a numerically safe way by
%  subtracting the row maximum to avoid cancelation after taking the exp
%  the sum is done along the rows.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.

function [y,x] = logsumexp2(logx)
  N = size(logx,2); max_logx = max(logx,[],2);
  % we have all values in the log domain, and want to calculate a sum
  x = exp(logx-max_logx*ones(1,N));
  y = log(sum(x,2)) + max_logx;