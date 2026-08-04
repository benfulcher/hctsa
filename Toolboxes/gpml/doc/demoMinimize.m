clear all, close all, nit = 50;     % tabular rasa, max number of fun/grad evals
n = 300; x = 10*sort(rand(n,1));                                   % sample data
sn = 0.2; y = sin(x)+sqrt(x) + randn(n,1)*sn;
cov = {@covSEiso}; sf = 1; ell = 0.4; hyp0.cov  = log([ell;sf]);    % covariance
mean = {@meanSum,{@meanLinear,@meanConst}}; hyp0.mean = [0.2;1];          % mean
lik = 'likGauss'; hyp0.lik = log(sn); inf = 'infExact';   % likelihood/inference
par = {inf, mean, cov, lik, x, y};


[x1 fx1 c1]   = minimize(        hyp0,@gp,-nit,par{:});
[x2 fx2 c2]   = minimize_minfunc(hyp0,@gp,-nit,par{:});
try
  [x3 fx3 c3] = minimize_lbfgsb(hyp0,@gp,-nit,par{:});
catch
  x3 = hyp0; fx3 = gp(hyp0,par{:}); c3 = 0;
end
fprintf('nit=%02d: minimize:          %1.2f\n',c1,fx1(end))
fprintf('nit=%02d: minimize_minfunc:  %1.2f\n',c2,fx2(end))
fprintf('nit=%02d: minimize_lbfgsb:   %1.2f\n',c3,fx3(end))