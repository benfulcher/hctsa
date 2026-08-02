disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
clear all, close all, write_fig = 0; N = 20;
sd = 3; rand('seed',sd), randn('seed',sd)       % set a seed for reproducability

fprintf('a) switch between FITC/VFE/SPEP via the opt.s parameter\n')
a = 0.3; b = 1.2; f = @(x) a*x + b + sin(x);               % underlying function
n = 30; sn = 0.5;          % number of training points, noise standard deviation
x = 2*rand(n,1)-1; x = 1+4*x+sign(x); y = f(x)+sn*randn(n,1);      % sample data
liktyp = input('Which likelihood?\n  (g)aussian, (l)ogistic: ','s');

cov = {@covSEiso}; sf = 2; ell = 1.0; hyp.cov = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; hyp.mean = [a;b];
if isequal(liktyp,'g')
  lik = {@likGauss};    hyp.lik = log(sn); inf = @infGaussLik;
else
  lik = {@likLogistic}; hyp.lik = [];      inf = @infLaplace;  y = sign(y);
end

fprintf('Optimise hyperparameters.\n')
hyp = minimize(hyp,@gp,-N,inf,mean,cov,lik,x,y);      % optimise hyperparameters
xs = linspace(-8,10,2e3)'; ys = f(xs);                   % exact function values

[ymu,ys2] = gp(hyp,inf,mean,cov,lik,x,y,xs);                  % dense prediction
[nlZ,dnlZ] = gp(hyp,inf,mean,cov,lik,x,y); % marginal likelihood and derivatives

% turn into a sparse approximation
xu = linspace(-3,5,25)'; cov = {'apxSparse',cov,xu};           % inducing points
infv  = @(varargin) inf(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
[ymuv,ys2v] = gp(hyp,infv,mean,cov,lik,x,y,xs);
infs = @(varargin) inf(varargin{:},struct('s',0.7));           % SPEP, 0<opt.s<1
[ymus,ys2s] = gp(hyp,infs,mean,cov,lik,x,y,xs);
inff = @(varargin) inf(varargin{:},struct('s',1.0));           % FITC, opt.s = 1
[ymuf,ys2f] = gp(hyp,inff,mean,cov,lik,x,y,xs);

fprintf('b) we can run sparse EP for FITC, as well\n')
infe = @infFITC_EP; [ymue,ys2e] = gp(hyp,infe,mean,cov,lik,x,y,xs);

subplot(211)
plot(xs,ymu,'k','LineWidth',2), hold on
plot(xs,ymuv,'g-.','LineWidth',2)
plot(xs,ymus,'m:','LineWidth',2)
plot(xs,ymuf,'c--','LineWidth',2)
plot(xs,ymue,'y:','LineWidth',2)
legend('exact','VFE','SPEP','FITC','FITC_E_P'), title('Predictive mean')
plot(x,y,'r+'), plot(xs,ys,'r')
plot(xs,ymu+2*sqrt(ys2),'k'), plot(xs,ymu-2*sqrt(ys2),'k')
xlim([-8,10]), ylim([-3,6])

subplot(212)
plot(xs,sqrt(ys2),'k','LineWidth',2), hold on
plot(xs,sqrt(ys2v),'g-.','LineWidth',2)
plot(xs,sqrt(ys2s),'m:','LineWidth',2)
plot(xs,sqrt(ys2f),'c--','LineWidth',2)
plot(xs,sqrt(ys2e),'y:','LineWidth',2)
legend('exact','VFE','SPEP','FITC','FITC_E_P'), title('Predictive standard dev')
xlim([-8,10]), if write_fig, print -depsc f11.eps; end

fprintf('c) specify inducing points via\n')
fprintf('1) hyp.xu or 2) {''apxSparse'',cov,xu}\n')
[nlZ1,dnlZ1] = gp(hyp,inf,mean,cov,lik,x,y); dnlZ1
hyp.xu = xu;
[nlZ2,dnlZ2] = gp(hyp,inf,mean,cov,lik,x,y); dnlZ2
fprintf('  The second has priority and\n')
fprintf('  results in derivatives w.r.t. xu\n')

fprintf('d) optimise nlZ w.r.t. inducing inputs\n')
fprintf('   by gradient descent\n')
hyp = minimize(hyp,@gp,-N,inf,mean,cov,lik,x,y);
if isequal(liktyp,'g')
  fprintf('   and by discrete swapping\n')
  z = 3*randn(100,1); % candidate inducing points
  nswap = N; % number of swaps between z and hyp.xu 
  [hyp,nlZ] = vfe_xu_opt(hyp,mean,cov,x,y,z,nswap);
end