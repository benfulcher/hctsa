% Optimize inducing inputs for the VFE approximation (not FITC).
%
% One can perform a gradient-based optimisation of the inducing inputs xu by
% specifying them via hyp.xu rather than through {@apxSparse,cov,xu}.
%
% An alternative way of optimising xu (in order to overcome local minima) is
% to simply compute the expected change in marginal likelihood of a set of
% candidate inducing points z and performing replacing the least relevant
% inducing input in xu with the most promising candidate from z. Efficient
% candidate scoring is only possible for the VFE approximation i.e. we use
% opt = struct('s',0.0); and lik = @likGauss; in the following.
%
% The call
%    [hyp,nlZ] = vfe_xu_opt(hyp,mean,cov,x,y, z,nswap);
% changes the inducing inputs in hyp.xu as to minimise nlZ. At most nswap
% swapping operations are performed between hyp.xu and the candidates in z.
% At the end the negative marginal likelihood nlZ is the same as obtained by
%     [post nlZ] = infGaussLik(hyp, mean, cov, lik, x, y, opt)
% where opt = struct('s',0.0); and lik = @likGauss.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2017-02-17.
%
% See also apx.m, infGaussLik.m.

function [hyp,nlZ] = vfe_xu_opt(hyp,mean,cov,x,y, z,nswap)

nlZmin = inf; u = hyp.xu; cov = cov{2};
for i=1:nswap
  [nlZ,du,dz] = vfe(hyp,mean,{'apxSparse',cov,u},x,y,z);
  if nlZ<nlZmin
    umin = u; nlZmin = nlZ;
    [dum,iu] = min(du); [dzm,iz] = min(dz);
    fprintf('%03d) nlZ=%1.4e -> %1.4e\n',i,nlZ,nlZ+dum+dzm)
    ui = u(iu,:); u(iu,:) = z(iz,:); z(iz,:) = ui;
  else
    fprintf('%03d) nlZ=%1.4e ~< %1.4e\n',i,nlZmin,nlZ)
    nlZ = nlZmin; u = umin; break
  end
end
hyp.xu = u;

function [nlZ,du,dz] = vfe(hyp,mean,cov,x,y,z)
  [n, D] = size(x);                                                 % dimensions
  m = feval(mean{:}, hyp.mean, x);                        % evaluate mean vector
  sn2 = exp(2*hyp.lik);                             % noise variance of likGauss
  xu = cov{3}; nu = size(xu,1);                        % extract inducing points
  k = @(x,z) feval(cov{2}{:},hyp.cov,x,z);    % shortcut for covariance function
  Kuu = k(xu,xu); Ku = k(xu,x); diagK = k(x,'diag');   % get the building blocks
  snu2 = 1e-6*(trace(Kuu)/nu);                 % stabilise by 0.1% of signal std
  Luu  = chol(Kuu+snu2*eye(nu));                       % Kuu + snu2*I = Luu'*Luu
  dgiKuu = sum((Luu\eye(nu)).^2,2);                             % diag(inv(Kuu))
  V  = Luu'\Ku; Lu = chol(eye(nu) + V*V'/sn2);    % V = inv(Luu')*Ku => V'*V = Q
  r = y-m; R = Luu\V; B = Lu'\V; alpha = r/sn2 - (B'*(B*r))/(sn2*sn2);
  t = max(diagK-sum(V.*V,1)',0)/sn2;
  nlZ = r'*alpha/2 + sum(log(diag(Lu))) + sum(t)/2 + n*log(2*pi*sn2)/2;
  if nargout>1
    q = R/sn2 - (R*B')*B/(sn2*sn2);
    du = (R*alpha).^2./(dgiKuu-sum(R.*q,2))/2   ... % aKa
                 + log(1-sum(R.*q,2)./dgiKuu)/2 ... % logdet Lui
                 + sum(R.*R,2)./dgiKuu/sn2/2;       % dst2
  end
  if nargout>2
    Kzu = k(z,xu); Kz = k(x,z); diagKz = k(z,'diag');    % query new cross-terms
    R = Kzu*R-Kz'; q = R/sn2 - (R*B')*B/(sn2*sn2);
    c = max(diagKz-sum((Kzu/Luu).^2,2),0);
    dz = -(R*alpha).^2./(c+sum(R.*q,2))/2 ...
                               + log(1+sum(R.*q,2)./c)/2 - sum(R.*R,2)./c/sn2/2;
  end