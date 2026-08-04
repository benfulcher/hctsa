function [K,dK] = covScale(cov, lsf, hyp, x, z)

% covScale - compose a covariance function as a scaled version of another
% one to model functions of the form f(x) = sf(x) f0(x), where sf(x) is a
% scaling function determining the function's standard deviation given f0(x)
% is normalised.
%
% The covariance function is parameterized as:
%     k(x,z) = sf(x) * k_0(x,z) * sf(z)
% with an important special case being
%     k(x,z) = sf^2  * k_0(x,z).
%
% You can either use K = covScale(cov, lsf, hyp, x, z) where the log scaling
% function lsf is a GPML mean function with hyperparameters hyp_sf yielding
%     hyp = [ hyp_cov
%             hyp_lsf ]
% as hyperparameters
% or you can use covScale(cov, hyp, x, z) to perform
% rescaling by a scalar value sf specified as an additional variable yielding
%     hyp = [ hyp_cov
%             log(sf) ]
% as hyperparameters.
%
% Copyright (c) by Carl Edward Rasmussen, Hannes Nickisch & Roman Garnett
%                                                                    2016-04-26.
%
% See also covFunctions.m.

if nargin==0, error('cov function must be specified'), end
if nargin<=1, lsf = []; end, narg = nargin;                % set a default value
if isnumeric(lsf)&&~isempty(lsf)  % shift parameters if sf contains actually hyp
  if nargin>3, z = x; end
  if nargin>2, x = hyp; end
  if nargin>1, hyp = lsf; end
  narg = nargin+1; lsf = [];
end

% below we us narg instead of nargin to be independent of the parameter shift
if isempty(lsf), ssf = '1'; else ssf = feval(lsf{:}); end     % number of hypers
if narg<4, K = ['(',feval(cov{:}),'+',ssf,')']; return, end
if narg<5, z = []; end                                     % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode
[n,D] = size(x);                                       % dimension of input data
ncov = eval(feval(cov{:})); hyp_cov = hyp(1:ncov); % number of params, split hyp
nsf  = eval(ssf);           hyp_lsf = hyp(ncov+(1:nsf));
scalar = isempty(lsf); if scalar, sf = exp(hyp_lsf); end
if ncov+nsf~=numel(hyp), error('Wrong number of hyper parameters.'), end

if nargout > 1
  [K0,dK0] = feval(cov{:},hyp_cov,x,z);
else
  K0 = feval(cov{:},hyp_cov,x,z);
end

if scalar
  sfx = sf; dsfx = @(q) sf*sum(q);
else
  [lsfx,dlsfx] = feval(lsf{:},hyp_lsf,x);
  sfx = exp(lsfx); dsfx = @(q) dlsfx(q.*sfx);
end
if dg
  S = sfx.*sfx; sfz = sfx; dsfz = dsfx;
else
  if xeqz
    sfz = sfx; dsfz = dsfx;
  else
    if scalar
      sfz = sf; dsfz = @(q) sf*sum(q);
    else
      [lsfz,dlsfz] = feval(lsf{:},hyp_lsf,z);
      sfz = exp(lsfz); dsfz = @(q) dlsfz(q.*sfz);
    end
  end
  S = sfx*sfz';
end

K = S.*K0;                                                          % covariance
if nargout > 1                                    % directional hyper derivative
  dK = @(Q) dirder(Q,S,K0,dK0,sfx,dsfx,sfz,dsfz,dg);
end

function [dhyp,dx] = dirder(Q,S,K0,dK0,sfx,dsfx,sfz,dsfz,dg)
  if nargout>1
    [dhyp0,dx0] = dK0(Q.*S); dx = dx0;            % sx contributions are missing
  else
    dhyp0 = dK0(Q.*S);
  end
  Q = Q.*K0;
  if dg
    qz = Q.*sfx; qx = Q.*sfz;
  else
    qz = sum(Q'*sfx,2); qx = sum(Q*sfz,2);
  end
  dhyp = [dhyp0; dsfx(qx)+dsfz(qz)];