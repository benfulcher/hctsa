function [K,dK] = covProd(cov, hyp, x, z)

% covProd - compose a covariance function as the product of other covariance
% functions. This function doesn't actually compute very much on its own, it
% merely does some bookkeeping, and calls other covariance functions to do the
% actual work.
%
% Note that cov = {cov1, cov2, .., false} turns of covariance matrix caching.
% This option slows down the computations but can help out if you products of
% many huge matrices lead to working memory shortage.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-18.
%
% See also covFunctions.m.

if isempty(cov), error('We require at least one factor.'), end
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
[n,D] = size(x); nh = numel(hyp);                    % number of hyperparameters

v = [];               % v vector indicates to which covariance parameters belong
for ii = 1:nc, v = [v repmat(ii, 1, eval(char(j(ii))))]; end

K = 1; Ki = cell(nc,1); dKi = cell(nc,1);                               % init K
for ii = 1:nc                                  % iteration over factor functions
  f = cov(ii); if iscell(f{:}), f = f{:}; end   % expand cell array if necessary
  if cache
    if nargin>1
      [Ki{ii},dKi{ii}] = feval(f{:}, hyp(v==ii), x, z); Kii = Ki{ii};    % track
    else
      Ki{ii} = feval(f{:}, hyp(v==ii), x, z); Kii = Ki{ii};
    end
  else
    Kii = feval(f{:}, hyp(v==ii), x, z);
  end
  K = K .* Kii;                                         % accumulate covariances
end
dK = @(Q) dirder(Q,Ki,dKi,v,nc,nh,cov,hyp,x,z,cache);   % directional derivative

function [dhyp,dx] = dirder(Q,Ki,dKi,v,nc,nh,cov,hyp,x,z,cache)
  dhyp = zeros(nh,1); dx = 0;
  for ii = 1:nc
    Qi = Q;
    for jj=1:nc     % accumulate
      if ii~=jj
        if cache
          Qi = Qi .* Ki{jj};
        else
          f = cov(jj); if iscell(f{:}), f = f{:}; end % expand cell if necessary
          Qi = Qi .* feval(f{:}, hyp(v==jj), x, z);
        end
      end
    end
    if cache
      dKii = dKi{ii};
    else
      f = cov(ii); if iscell(f{:}), f = f{:}; end     % expand cell if necessary
      [junk,dKii] = feval(f{:}, hyp(v==ii), x, z);
    end
    if nargout==1
      dhyp(v==ii,1) = dKii(Qi);
    else
      [dhyp(v==ii,1),dxi] = dKii(Qi); dx = dx+dxi;
    end
  end