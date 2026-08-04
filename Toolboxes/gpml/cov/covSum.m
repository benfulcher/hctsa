function [K,dK] = covSum(cov, hyp, x, z)

% covSum - compose a covariance function as the sum of other covariance
% functions. This function doesn't actually compute very much on its own, it
% merely does some bookkeeping, and calls other covariance functions to do the
% actual work.
%
% Note that cov = {cov1, cov2, .., false} turns of covariance matrix caching.
% This option slows down the computations but can help out if you sums of
% many huge matrices lead to working memory shortage.
%
% Copyright (c) by Carl Edward Rasmussen & Hannes Nickisch 2016-04-18.
%
% See also covFunctions.m.

if isempty(cov), error('We require at least one summand.'), end
if isnumeric(cov{end}) || islogical(cov{end})          % detect whether to cache
  cache = ~isequal(0,cov{end}) && ~isequal(false,cov{end});
  cov = cov(1:end-1);                           % chop off last element from cov
else
  cache = true;
end
nc = numel(cov);                        % number of terms in covariance function
for ii = 1:nc                                % iterate over covariance functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  j(ii) = cellstr(feval(f{:}));                          % collect number hypers
end

if nargin<3                                        % report number of parameters
  K = char(j(1)); for ii=2:nc, K = [K, '+', char(j(ii))]; end, return
end
if nargin<4, z = []; end                                   % make sure, z exists
[n,D] = size(x);

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:nc, v = [v repmat(ii, 1, eval(char(j(ii))))]; end

K = 0; dKi = cell(nc,1);                                                % init K
for ii = 1:nc                                 % iteration over summand functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  if nargout>1 && cache
    [Kii,dKi{ii}] = feval(f{:}, hyp(v==ii), x, z);                  % keep track
  else
    Kii = feval(f{:}, hyp(v==ii), x, z);
  end
  K = K + Kii;                                          % accumulate covariances
end
dK = @(Q) dirder(Q,dKi,v,nc,cov,hyp,x,z,cache);         % directional derivative

function [dhyp,dx] = dirder(Q,dKi,v,nc,cov,hyp,x,z,cache)
  dhyp = zeros(size(v,2),1); dx = 0;
  for ii = 1:nc
    if cache
      dKii = dKi{ii};
    else
      f = cov(ii); if iscell(f{:}), f = f{:}; end     % expand cell if necessary
      [junk,dKii] = feval(f{:}, hyp(v==ii), x, z);
    end
    if nargout > 1
      [dhyp(v==ii,1),dxi] = dKii(Q); dx = dx+dxi;
    else
      dhyp(v==ii,1) = dKii(Q);
    end
  end