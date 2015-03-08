function K = covScale(cov, lsf, hyp, x, z, i)

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
% You can either use K = covScale(cov, lsf, hyp, x, z, i) where the log scaling
% function lsf is a GPML mean function with hyperparameters hyp_sf yielding
%     hyp = [ hyp_lsf
%             hyp_cov ]
% as hyperparameters
% or you can use covScale(cov, hyp, x, z, i) to perform
% rescaling by a scalar value sf specified as an additional variable yielding
%     hyp = [ log(sf)
%             hyp_cov ]
% as hyperparameters.
%
% Copyright (c) by Carl Edward Rasmussen, Hannes Nickisch & Roman Garnett
%                                                                    2014-09-05.
%
% See also COVFUNCTIONS.M.

if nargin==0, error('cov function must be specified'), end
if nargin<=1, lsf = []; end, narg = nargin;                % set a default value
if isnumeric(lsf)&&~isempty(lsf)  % shift parameters if sf contains actually hyp
  if nargin>4, i = z; end
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
nsf  = eval(ssf);           hyp_lsf = hyp(1:nsf);  % number of params, split hyp
ncov = eval(feval(cov{:})); hyp_cov = hyp(nsf+(1:ncov));
scalar = isempty(lsf); if scalar, sf = exp(hyp_lsf); end

if scalar, sfx = sf; else sfx = exp(feval(lsf{:},hyp_lsf,x)); end
if dg
  K = sfx.*sfx;
else
  if xeqz, sfz = sfx;
  else if scalar, sfz = sf; else sfz = exp(feval(lsf{:},hyp_lsf,z)); end, end
  K = sfx*sfz';
end

if narg<6                                                          % covariances
  K = K.*feval(cov{:},hyp_cov,x,z);
else                                                               % derivatives
  if i>nsf                             % wrt covariance function hyperparameters
    K = K.*feval(cov{:},hyp_cov,x,z,i-nsf);
  else                                    % wrt scaling function hyperparameters
    K = feval(cov{:},hyp_cov,x,z);
    if scalar, dsfx = sfx; else dsfx = sfx.*feval(lsf{:},hyp_lsf,x,i); end
    if dg
      K = (2*sfx.*dsfx).*K;
    else
      if xeqz, dsfz = dsfx;
      else
        if scalar, dsfz = sfz; else dsfz = sfz.*feval(lsf{:},hyp_lsf,z,i); end
      end
      K = (sfx*dsfz'+dsfx*sfz').*K;
    end
  end
end