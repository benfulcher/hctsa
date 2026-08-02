function [m,dm] = meanGPexact(mean,cov,x,y, hypz,z)

% Mean function being the predictive mean of a GP model:
%
% mu(z) = posterior mean of GP at location z as given by
% mu(z) = gp(hyp,@infExact,mean,cov,@likGauss,x,y, z) where
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
% Copyright (c) by Hannes Nickisch, 2016-04-16.
%
% See also meanFunctions.m and meanGP.m.

if nargin<4, error('GP must be specified.'), end           % check for dimension
if isempty(mean), mean = @meanZero; end              % set default and make cell
if ~iscell(mean), mean = {mean};    end
if isempty(cov),  cov  = @covSEiso; end
if ~iscell(cov),  cov  = {cov};     end
nms = feval(mean{:}); ncs = feval(cov{:});     % number of hyperparameter string
if nargin<6, m = [ncs,'+1+',nms]; return, end % report number of hyperparameters

[nz,D] = size(z); n = size(x,1);
nc = eval(ncs); nm = eval(nms);
hyp = vec2any(struct('cov',zeros(nc,1),'lik',0,'mean',zeros(nm,1)),hypz);

[mu,dmu]  = feval(mean{:},hyp.mean,x);
[muz,dmuz] = feval(mean{:},hyp.mean,z);
[K,dK]  = feval(cov{:}, hyp.cov, x);
[kz,dkz] = feval(cov{:}, hyp.cov, x,z);

sn2 = exp(2*hyp.lik);                               % noise variance of likGauss
if sn2<1e-6                        % very tiny sn2 can lead to numerical trouble
  L = chol(K+sn2*eye(n)); sl =   1;   % Cholesky factor of covariance with noise
else
  L = chol(K/sn2+eye(n)); sl = sn2;                       % Cholesky factor of B
end
iKs = @(t) solve_chol(L,t)/sl;                       % iKs(t) = (K+sn2*eye(n))\t
alpha = iKs(y-mu);
m = muz+kz'*alpha;                                         % eval posterior mean
dm = @(q) dirder(q,alpha,dmu,dmuz,kz,dkz,iKs,dK,sn2);   % directional derivative

function dmdhyp = dirder(q,alpha,dmu,dmuz,kz,dkz,iKs,dK,sn2)
  v = iKs(kz*q);
  dmdhyp = [dkz(alpha*q')-dK(alpha*v'); -2*sn2*v'*alpha; dmuz(q)-dmu(v)];