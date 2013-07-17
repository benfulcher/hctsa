function A = meanSum(mean, hyp, x, i)

% meanSum - compose a mean function as the sum of other mean functions.
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean functions to do the actual work.
%
% m(x) = \sum_i m_i(x)
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2010-08-04.
%
% See also MEANFUNCTIONS.M.

for ii = 1:numel(mean)                             % iterate over mean functions
  f = mean(ii); if iscell(f{:}), f = f{:}; end  % expand cell array if necessary
  j(ii) = cellstr(feval(f{:}));                          % collect number hypers
end

if nargin<3                                        % report number of parameters
  A = char(j(1)); for ii=2:length(mean), A = [A, '+', char(j(ii))]; end; return
end

[n,D] = size(x);

v = [];                     % v vector indicates to which mean parameters belong
for ii = 1:length(mean), v = [v repmat(ii, 1, eval(char(j(ii))))]; end

if nargin==3                                               % compute mean vector
  A = zeros(n,1);                                               % allocate space
  for ii = 1:length(mean)                     % iteration over summand functions
    f = mean(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if needed
    A = A + feval(f{:}, hyp(v==ii), x);                       % accumulate means
  end
else                                                 % compute derivative vector
  if i<=length(v)
    ii = v(i);                                             % which mean function
    j = sum(v(1:i)==ii);                          % which parameter in that mean
    f = mean(ii);
    if iscell(f{:}), f = f{:}; end         % dereference cell array if necessary
    A = feval(f{:}, hyp(v==ii), x, j);                      % compute derivative
  else
    A = zeros(n,1);
  end
end