function [m,dm] = meanProd(mean, hyp, x)

% meanProd - compose a mean function as the product of other mean functions.
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls other mean functions to do the actual work.
%
% m(x) = \prod_i m_i(x)
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

m = ones(n,1); mi = cell(nm,1); dmi = cell(nm,1);               % allocate space
for ii = 1:nm                                  % iteration over factor functions
  f = mean(ii); if iscell(f{:}), f = f{:}; end     % expand cell array if needed
  [mi{ii},dmi{ii}] = feval(f{:}, hyp(v==ii), x);
  m = m.*mi{ii};                                              % accumulate means
end
dm = @(q) dirder(q,mi,dmi,v,nm);                        % directional derivative

function dhyp = dirder(q,mi,dmi,v,nm)
  dhyp = zeros(size(v,2),1);
  for ii = 1:nm
    qi = q; for jj=1:nm, if ii~=jj, qi = qi .* mi{jj}; end, end     % accumulate
    dhyp(v==ii,1) = dmi{ii}(qi);
  end