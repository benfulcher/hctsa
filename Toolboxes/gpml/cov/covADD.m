function [K,dK] = covADD(cov, hyp, x, z)

% Additive covariance function using 1d base covariance functions 
% cov_d(x_d,z_d;hyp_d) with individual hyperparameters hyp_d, d=1..D.
%
% k  (x,z) = \sum_{r \in R} sf^2_r k_r(x,z), where 1<=r<=D and
% k_r(x,z) = \sum_{|I|=r} \prod_{i \in I} cov_i(x_i,z_i;hyp_i)
%
% hyp = [ hyp_1
%         hyp_2
%          ...
%         hyp_D 
%         log(sf_R(1))
%          ...
%         log(sf_R(end)) ]
%
% where hyp_d are the parameters of the 1d covariance function which are shared
% over the different values of r=R(1),..,R(end) where 1<=r<=D.
%
% Usage: covADD({[1,3,4],cov}, ...) or
%        covADD({[1,3,4],cov_1,..,cov_D}, ...).
%
% Please see the paper "Additive Gaussian Processes" by Duvenaud, Nickisch and 
% Rasmussen, NIPS, 2011 for details.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2017-02-21.
%                  multiple covariance support contributed by Truong X. Nghiem
%
% See also covFunctions.m.

R = fix(cov{1});                            % only positive integers are allowed
if min(R)<1, error('only positive R up to D allowed'), end
nr = numel(R);                      % number of different degrees of interaction
nc = numel(cov)-1; D = 1;         % number of provided covariances, D for indiv.
for j=1:nc, if ~iscell(cov{j+1}), cov{j+1} = {cov{j+1}}; end, end
if nc==1, nh = eval(feval(cov{2}{:}));  % no of hypers per individual covariance
else nh = zeros(nc,1); for j=1:nc, nh(j) = eval(feval(cov{j+1}{:})); end, end

if nargin<3                                  % report number of hyper parameters
  if nc==1, K = ['D*', int2str(nh), '+', int2str(nr)];
  else
    K = ['(',int2str(nh(1))]; for ii=2:nc, K = [K,'+',int2str(nh(ii))]; end
    K = [K, ')+', int2str(nr)];
  end
  return
end
if nargin<4, z = []; end                                   % make sure, z exists

[n,D] = size(x);                                                % dimensionality
if nc==1, nh = ones(D,1)*nh; cov = [cov(1),repmat(cov(2),1,D)]; end
nch = sum(nh);                                      % total number of cov hypers
sf2 = exp( 2*hyp(nch+(1:nr)) );         % signal variances of individual degrees

[Kd,dKd] = Kdim(cov(2:end),nh,hyp(1:nch),x,z);  % eval dimensionwise covariances
EE = elsympol(Kd,max(R));                 % Rth elementary symmetric polynomials
K = 0; for ii=1:nr, K = K + sf2(ii)*EE(:,:,R(ii)+1); end      % sf2 weighted sum

if nargout > 1
  dK = @(Q) dirder(Q,Kd,dKd,EE,R,sf2);
end

function [dhyp,dx] = dirder(Q,Kd,dKd,EE,R,sf2)
  D = numel(dKd); n = size(Q,1); nr = numel(R); dhyp = zeros(0,1);
  if nargout > 1, dx = zeros(n,D); end                    % allocate if required
  for j=1:D
    % the final derivative is a sum of multilinear terms, so if only one term
    % depends on the hyperparameter under consideration, we can factorise it 
    % out and compute the sum with one degree less, the j-th elementary
    % covariance depends on the hyperparameter batch hyp(j) and inputs x(:,j)
    E = elsympol(Kd(:,:,[1:j-1,j+1:D]),max(R)-1);  %  R-1th elementary sym polyn
    Qj = 0; for ii=1:nr, Qj = Qj + sf2(ii)*E(:,:,R(ii)); end  % sf2 weighted sum
    if nargout > 1, [dhypj,dxj] = dKd{j}(Qj.*Q); dx(:,j) = dxj;
    else dhypj = dKd{j}(Qj.*Q); end, dhyp = [dhyp; dhypj];
  end
  dhyp = [dhyp; 2*squeeze(sum(sum(bsxfun(@times,EE(:,:,R+1),Q),1),2)).*sf2];

% evaluate dimensionwise covariances K
function [K,dK] = Kdim(cov,nh,hyp,x,z)
  [n,D] = size(x);                                              % dimensionality
  xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;      % determine mode

  if dg                                                        % allocate memory
    K = zeros(n,1,D);
  else
    if xeqz, K = zeros(n,n,D); else K = zeros(n,size(z,1),D); end
  end
  dK = cell(D,1);

  for d=1:D
    hyp_d = hyp(sum(nh(1:d-1))+(1:nh(d)));        % hyperparamter of dimension d
    if dg
      [K(:,:,d),dK{d}] = feval(cov{d}{:},hyp_d,x(:,d),'diag');
    else
      if xeqz
        [K(:,:,d),dK{d}] = feval(cov{d}{:},hyp_d,x(:,d));
      else
        [K(:,:,d),dK{d}] = feval(cov{d}{:},hyp_d,x(:,d),z(:,d));
      end
    end
  end