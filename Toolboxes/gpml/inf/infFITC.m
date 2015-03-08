function [post nlZ dnlZ] = infFITC(hyp, mean, cov, lik, x, y)

% FITC approximation to the posterior Gaussian process. The function is
% equivalent to infExact with the covariance function:
%   Kt = Q + G; G = diag(g); g = diag(K-Q);  Q = Ku'*inv(Quu)*Ku;
% where Ku and Kuu are covariances w.r.t. to inducing inputs xu, snu2 = sn2/1e6
% is the noise of the inducing inputs and Quu = Kuu + snu2*eye(nu).
% We fixed the standard deviation of the inducing inputs snu to be a one per mil
% of the measurement noise's standard deviation sn.
% The implementation exploits the Woodbury matrix identity
%   inv(Kt) = inv(G) - inv(G)*V'*inv(eye(nu)+V*inv(G)*V')*V*inv(G)
% in order to be applicable to large datasets. The computational complexity
% is O(n nu^2) where n is the number of data points x and nu the number of
% inducing inputs in xu.
%
% The function takes a specified covariance function (see covFunctions.m) and
% likelihood function (see likFunctions.m), and is designed to be used with
% gp.m and in conjunction with covFITC and likGauss.
%
% The inducing points can be specified through 1) the 2nd covFITC parameter or
% by 2) providing a hyp.xu hyperparameters. Note that 2) has priority over 1).
% In case 2) is provided and derivatives dnlZ are requested, there will also be
% a dnlZ.xu field allowing to optimise w.r.t. to the inducing points xu. However
% the derivatives dnlZ.xu can only be computed for one of the following eight
% covariance functions: cov{Matern|PP|RQ|SE}{iso|ard}.
%
% Copyright (c) by Ed Snelson, Carl Edward Rasmussen 
%                                               and Hannes Nickisch, 2013-10-28.
%
% See also INFMETHODS.M, COVFITC.M.

if iscell(lik), likstr = lik{1}; else likstr = lik; end
if ~ischar(likstr), likstr = func2str(likstr); end
if ~strcmp(likstr,'likGauss')               % NOTE: no explicit call to likGauss
  error('Inference with inFITC only possible with Gaussian likelihood.');
end
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if ~strcmp(cov1,'covFITC'); error('Only covFITC supported.'), end    % check cov
if isfield(hyp,'xu'), cov{3} = hyp.xu; end  % hyp.xu is provided, replace cov{3}

[diagK,Kuu,Ku] = feval(cov{:}, hyp.cov, x);         % evaluate covariance matrix
m = feval(mean{:}, hyp.mean, x);                          % evaluate mean vector
[n, D] = size(x); nu = size(Kuu,1);
cov{3} = reshape(cov{3},[],D);

sn2  = exp(2*hyp.lik);                              % noise variance of likGauss
snu2 = 1e-6*sn2;                              % hard coded inducing inputs noise
Luu  = chol(Kuu+snu2*eye(nu));                         % Kuu + snu2*I = Luu'*Luu
V  = Luu'\Ku;                                     % V = inv(Luu')*Ku => V'*V = Q
g_sn2 = diagK + sn2 - sum(V.*V,1)';          % g + sn2 = diag(K) + sn2 - diag(Q)
Lu = chol(eye(nu) + (V./repmat(g_sn2',nu,1))*V');  % Lu'*Lu=I+V*diag(1/g_sn2)*V'
r  = (y-m)./sqrt(g_sn2);
be = Lu'\(V*(r./sqrt(g_sn2)));
iKuu = solve_chol(Luu,eye(nu));                       % inv(Kuu + snu2*I) = iKuu
post.alpha = Luu\(Lu\be);                      % return the posterior parameters
post.L  = solve_chol(Lu*Luu,eye(nu)) - iKuu;                    % Sigma-inv(Kuu)
post.sW = ones(n,1)/sqrt(sn2);           % unused for FITC prediction  with gp.m

if nargout>1                                % do we want the marginal likelihood
  nlZ = sum(log(diag(Lu))) + (sum(log(g_sn2)) + n*log(2*pi) + r'*r - be'*be)/2;
  if nargout>2                                         % do we want derivatives?
    dnlZ = hyp;                                 % allocate space for derivatives
    al = r./sqrt(g_sn2) - (V'*(Lu\be))./g_sn2;      % al = (Kt+sn2*eye(n))\(y-m)
    B = iKuu*Ku; w = B*al;
    W = Lu'\(V./repmat(g_sn2',nu,1));
    for i = 1:numel(hyp.cov)
      [ddiagKi,dKuui,dKui] = feval(cov{:}, hyp.cov, x, [], i);  % eval cov deriv
      R = 2*dKui-dKuui*B; v = ddiagKi - sum(R.*B,1)';   % diag part of cov deriv
      dnlZ.cov(i) = (ddiagKi'*(1./g_sn2)+w'*(dKuui*w-2*(dKui*al))-al'*(v.*al)...
                                 - sum(W.*W,1)*v - sum(sum((R*W').*(B*W'))) )/2;
    end
    diag_dK = 1./g_sn2 - sum(W.*W,1)' - al.*al;                  % diag(dnlZ/dK)
    dnlZ.lik = sn2*sum(diag_dK);
    % since snu2 is a fixed fraction of sn2, there is a covariance-like term in
    BWt = B*W'; % the derivative as well
    dKuui = 2*snu2; R = -dKuui*B; v = -sum(R.*B,1)';    % diag part of cov deriv
    dnlZ.lik = dnlZ.lik + (w'*dKuui*w -al'*(v.*al)...
                                    - sum(W.*W,1)*v - sum(sum((R*W').*BWt)) )/2; 
    for i = 1:numel(hyp.mean)
      dnlZ.mean(i) = -feval(mean{:}, hyp.mean, x, i)'*al;
    end
    if isfield(hyp,'xu')                 % derivatives w.r.t. inducing points xu
      xu = cov{3};
      cov = cov{2};           % get the non FITC part of the covariance function
      Kpu  = cov_deriv_sq_dist(cov,hyp.cov,xu,x);           % d K(xu,x ) / d D^2
      Kpuu = cov_deriv_sq_dist(cov,hyp.cov,xu);             % d K(xu,xu) / d D^2
      if iscell(cov), covstr = cov{1}; else covstr = cov; end
      if ~ischar(covstr), covstr = func2str(covstr); end
      if numel(strfind(covstr,'iso'))>0            % characteristic length scale
        e = 2*exp(-2*hyp.cov(1));
      else
        e = 2*exp(-2*hyp.cov(1:D));
      end
      v = diag_dK-1./g_sn2;        % BdK = B * ( dnlZ/dK - diag(diag(dnlZ/dK)) )
      BdK = B.*repmat(v',nu,1) + BWt*W + (B*al)*al';
      A = Kpu.*BdK; C = Kpuu.*(BdK*B'); C = diag(sum(C,2)-sum(A,2)) - C;
      dnlZ.xu = A*x*diag(e) + C*xu*diag(e);  % bring in data and inducing points
    end
  end
end
