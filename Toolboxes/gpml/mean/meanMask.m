function [m,dm] = meanMask(mask, mean, hyp, x)
% MEANMASK Apply a mean function to a subset of the dimensions only.
%  [m, dm] = MEANMASK (mask, mean, hyp, x)
%
% The subset of dimensions is specified by the input argument mask, and can
% either be a 0/1 mask, by a boolean mask, or by an index set. The second input
% argument, mean, should be a handle to a valid gpml mean function. The other
% input arguments are the same as for regular mean functions. The size of
% the hyperparamter vector, hyp, should match the number of hyperparameters of
% the masked mean function.
%
% For example, if the input x is 3-dimensional and we want a linear mean on the
% second component of the inputs we could write any of the following
%
% linear_on_x2 = {@meanMask, [0,1,0], @meanLinear}
% linear_on_x2 = {@meanMask, [false,true,false], @meanLinear}
% linear_on_x2 = {@meanMask, 2, @meanLinear}
%
% the function linear_on_x2 behaves like a regular mean function but the number
% of hyperparameters is 1, instead of 3.
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls another mean function to do the actual work.
%
% See also MEANFUNCTIONS

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-16.

if ~iscell(mean), mean = {mean}; end             % properly wrap into cell array
nh_string = feval(mean{:});         % number of hyperparameters of the full mean

if max(mask)<2 && length(mask)>1, mask = find(mask); end    % convert 1/0->index
D = length(mask);                                             % masked dimension
if nargin<4, m = num2str(eval(nh_string)); return; end    % number of parameters

if eval(nh_string)~=length(hyp)                          % check hyperparameters
  error('number of hyperparameters does not match size of masked data')
end

[m,dm] = feval(mean{:}, hyp, x(:,mask));
