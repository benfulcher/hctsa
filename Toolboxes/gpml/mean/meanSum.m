function [m,dm] = meanSum(mean, hyp, x)

% meanSum - compose a mean function as the sum of other mean functions.
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean functions to do the actual work.
%
% m(x) = \sum_i m_i(x)
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2020-07-14.
%
% See also meanFunctions.m.

nm = numel(mean);
for ii = 1:nm                                      % iterate over mean functions
  f = mean(ii); if iscell(f{:}), f = f{:}; end  % expand cell array if necessary
  j(ii) = cellstr(feval(f{:}));                          % collect number hypers
end

if nargin<3                                        % report number of parameters
  m = char(j(1)); for ii=2:nm, m = [m, '+', char(j(ii))]; end; return
end
[n,D] = size(x);

v = [];                     % v vector indicates to which mean parameters belong
for ii = 1:nm, v = [v repmat(ii, 1, eval(char(j(ii))))]; end

m = zeros(n,1); dmi = cell(nm,1);                               % allocate space
for ii = 1:nm                                 % iteration over summand functions
  f = mean(ii); if iscell(f{:}), f = f{:}; end     % expand cell array if needed
  [mi,dmi{ii}] = feval(f{:}, hyp(v==ii), x);
  m = m+mi;                                                   % accumulate means
end
dm = @(q) dirder(q,dmi,v,nm);                           % directional derivative

function dhyp = dirder(q,dmi,v,nm)
  dhyp = zeros(size(v,2),1);
  for ii = 1:nm, dhyp(v==ii,1) = dmi{ii}(q); end