clear all, close all

% GP spec
mean = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
cov = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
lik = @likGauss; sn = 0.1; hyp.lik = log(sn);
inf = @infGaussLik;

% synthetic data
n = 20; x = gpml_randn(0.3, n, 1);
K = feval(cov{:}, hyp.cov, x);
m = feval(mean{:}, hyp.mean, x);
y = chol(K)'*gpml_randn(0.15, n, 1) + m + exp(hyp.lik)*gpml_randn(0.2, n, 1);

% inference
[nlZ,dnlZ,post] = gp(hyp, inf, mean, cov, lik, x, y);

% prediction
ns = 101; xs = linspace(-1.9, 1.9, ns)';
[ymu,ys2,fmu,fs2] = gp(hyp, inf, mean, cov, lik, x, post, xs);

% compute cross-covariances
Lchol = isnumeric(post.L) && all(all(tril(post.L,-1)==0)&diag(post.L)'>0&isreal(diag(post.L))');
if Lchol
  Kss = feval(cov{:}, hyp.cov, xs);
  Ks = feval(cov{:}, hyp.cov, x, xs);
  V = post.L'\(repmat(post.sW,1,ns).*Ks);
  S2 = Kss - V'*V;

  assert( norm(diag(S2)-fs2)<1e-10 )
end