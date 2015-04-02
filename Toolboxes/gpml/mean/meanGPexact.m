function A = meanGPexact(mean,cov,x,y, hypz,z,i)

% Mean function being the predictive mean of a GP model:
%
% m(z) = posterior mean of GP at location z as given by
% m(z) = gp(hyp,@infExact,mean,cov,@likGauss,x,y, z) where
% hyp.mean = hyp_mean; hyp.lik = log(sn); hyp.cov = hyp.cov;
%
% The hyperparameters are:
%
% hypz = [ hyp_cov
%          log(sn)
%          hyp_mean ]
%
% where hyp_cov are the covariance function hyperparameters, sn is the
% noise variance of the Gaussian likelihood and hyp_mean are the mean
% function hyperparameters.
%
% Copyright (c) by Hannes Nickisch, 2014-11-01.
%
% See also MEANFUNCTIONS.M and MEANGP.M.

if nargin<4, error('GP must be specified.'), end           % check for dimension
if isempty(mean), mean = @meanZero; end              % set default and make cell
if ~iscell(mean), mean = {mean};    end
if isempty(cov),  cov  = @covSEiso; end
if ~iscell(cov),  cov  = {cov};     end
nms = feval(mean{:}); ncs = feval(cov{:});     % number of hyperparameter string
if nargin<6, A = [ncs,'+1+',nms]; return, end % report number of hyperparameters

[nz,D] = size(z); n = size(x,1);
nc = eval(ncs); nm = eval(nms);
hyp = rewrap(struct('cov',zeros(nc,1),'lik',0,'mean',zeros(nm,1)),hypz);

m  = feval(mean{:},hyp.mean,x);
mz = feval(mean{:},hyp.mean,z);
K  = feval(cov{:}, hyp.cov, x);
kz = feval(cov{:}, hyp.cov, x,z);

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
  L = chol(K+sn2*eye(n)); sl =   1;   % Cholesky factor of covariance with noise
else
  L = chol(K/sn2+eye(n)); sl = sn2;                       % Cholesky factor of B
end
iKs = @(t) solve_chol(L,t)/sl;                       % iKs(t) = (K+sn2*eye(n))\t
alpha = iKs(y-m);
if nargin==6                                               % eval posterior mean
  A = mz+kz'*alpha;
else
  if i<=nc                                          % covariance function hypers
    dK = feval(cov{:},hyp.cov,x,[],i);
    dkz = feval(cov{:},hyp.cov,x,z,i);
    A = dkz'*alpha-kz'*iKs(dK*alpha);
  elseif i==nc+1                              % likelihood function parameter sn
    A = -2*sn2*kz'*iKs(alpha);
  else                                                 % mean function parameter
    dm  = feval(mean{:},hyp.mean,x,i-nc-1);
    dmz = feval(mean{:},hyp.mean,z,i-nc-1);
    A = dmz-kz'*iKs(dm);
  end
end
