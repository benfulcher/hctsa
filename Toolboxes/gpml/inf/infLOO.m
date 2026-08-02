function [post nloo dnloo loo] = infLOO(hyp, mean, cov, lik, x, y)

% Leave-One-Out - Perform Least-Squares GP predictions in the style of the 
% Pseudo Likelihood in 5.4.2 p. 117 and the Probabilistic Least-squares 
% Classifier of 6.5  p. 146 of the GPML book.
% For Gaussian likelihood, the prediction is the same as with infExact, for
% lik{Erf,Logistic} the result corresponds to probabilistic least-squares
% classification. In any case, the objective returned is not the marginal
% likelihood but the predictive log probability.
% The standard deviation of the noise either comes from the sn parameter of the 
% likelihood function lik{Gauss,Laplace,Sech2} or from the covariance function
% by means of using covSum and covNoise i.e. for lik{Uni,Erf,Logistic}.
% We then compute the negative leave-one-out predictive probability and its 
% derivatives w.r.t. the hyperparameters. See also "help infMethods".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-05-02
%
% See also infMethods.m.

[n, D] = size(x);
if isnumeric(cov),  K = cov;                    % use provided covariance matrix
else [K,dK] = feval(cov{:},  hyp.cov,  x); end     % covariance matrix and deriv
if isnumeric(mean), m = mean;                         % use provided mean vector
else [m,dm] = feval(mean{:}, hyp.mean, x); end           % mean vector and deriv
if numel(hyp.lik)==0
  cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
  if ~strcmp(cov1,'covSum'); error('Only covSum supported.'), end    % check cov
  cov2 = cov{2}; npar = 0;
  for i=1:length(cov2)
    cov2i = cov2{i};
    if numel(cov2i)==1 || ischar(cov2i)
      if isa(cov2i, 'function_handle'), cov2i = func2str(cov2i); end
      if strcmp(cov2i,'covNoise'), isn = npar+1; sn2 = exp(2*hyp.cov(isn)); end
    end
    npar = npar + eval(feval(cov2i));
  end
  if ~exist('sn2','var'), error('We need covNoise in the covSum.'), end 
else                                                     % Gauss, Laplace, Sech2
  sn2 = exp(2*hyp.lik(end)); isn = 0;
  K   = K + sn2*eye(n);
end

L = chol(K)/sqrt(sn2);                % Cholesky factor of covariance with noise
alpha = solve_chol(L,y-m)/sn2;

post.alpha = alpha;                            % return the posterior parameters
post.sW = ones(n,1)/sqrt(sn2);                  % sqrt of noise precision vector
post.L  = L;                                        % L = chol(eye(n)+sW*sW'.*K)

if nargout>1                               % do we want the marginal likelihood?
  sdeff = 1./post.sW;                                 % effective Gaussian width
  mneff = K*alpha + m;     % alpha = inv(K)*(mneff-m), effective Gaussian center
  % all loo marginal means and variances at once
  V = post.L'\diag(1./sdeff); r = sum(V.*V,1)';               % r=diag( inv(K) )
  loo.fmu = mneff - alpha./r;                              % GPML book, eq. 5.12
  loo.fs2 = 1./r - sn2;                                    % GPML book, eq. 5.12
  [loo.lp, loo.ymu, loo.ys2] = feval(lik{:},hyp.lik,y,loo.fmu,loo.fs2);
  nloo = -sum(loo.lp);           % negative leave-one-out predictive probability
  if nargout>2                                         % do we want derivatives?
    [lZ,dlZmn,d2lZmn] = feval(lik{:},hyp.lik,y,loo.fmu,loo.fs2,'infEP');
    dlZva = (d2lZmn+dlZmn.*dlZmn)/2;                % derivative w.r.t. variance
    dnloo = hyp;
    iK = solve_chol(L,eye(n))/sn2;                % based on GPML book, eq. 5.13
    Q = iK*diag((dlZmn.*alpha-dlZva)./(r.*r))*iK - (iK*(dlZmn./r))*alpha';
    dnloo.cov = dK(Q);
    if isn, dnloo.cov(isn) = dnloo.cov(isn) + 2*sn2*sum(dlZva); end
    for j=1:numel(hyp.lik)     
      dnloo.lik(j) = -sum( feval(lik{:},hyp.lik,y,loo.fmu,loo.fs2,'infEP',j) );
      if j==numel(hyp.lik)
        dva = 2*sn2*sum(iK.*iK,2)./(r.*r);                 % GPML book, eq. 5.13
        dmn = 2*solve_chol(L,alpha)./r - alpha.*dva;       % GPML book, eq. 5.13
        dva = dva - 2*sn2;
        dnloo.lik(j) = dnloo.lik(j) -dlZmn'*dmn -dlZva'*dva;
      end
    end
    dnloo.mean = -dm(solve_chol(L,dlZmn./(sn2*r)));
  end
end
