clear all, close all                                                   % tidy up

n = 2350; x = 7*rand(n,1); xs = linspace(-0.4,8,200)'; opt = [];
xs(13:15) = xs(12);                                     % mess with the test set
x(14) = x(15); x(17) = x(16);                       % mess with the training set
sn = 0.05; ell = 2.7; sf = 1.3; y = 6*sin(pi*x)./(pi*x+1) + sn*randn(size(x));

hyp.lik = log(sn); lik = @likGauss; inf = @infGaussLik;             % likelihood
mean = {'meanSum',{'meanLinear','meanConst'}}; hyp.mean = [0.05;0.1];     % mean
cov = {'covMaterniso',7}; hyp.cov = log([ell;sf]);                  % covariance

mode = input('Which mode?\n (0) Regression, (1-4) Classification? ','s');
mode = str2double(mode);
switch mode
  case 1, lik = @likErf;      inf = @infLaplace;
  case 2, lik = @likErf;      inf = @infEP;
  case 3, lik = @likLogistic; inf = @infVB;
  case 4, lik = @likLogistic; inf = @infKL;
end
if mode>0, hyp.lik = []; y = sign(y); opt.inf = inf; end
fprintf('Likelihood is %s, inference method is %s.\n',func2str(lik),func2str(inf))

opt.ngrid = 200;     % number of equispaced grid points for matrix interpolation
% opt.deriv = false;                      % in case derivatives are not required
% opt.balance = true;     % balance state space for improved numerical stability
% opt.qridge = 1e-5; % add a ridge to the Q submatrices to improve the condition
% opt.rridge = 1e-5;                 % add a ridge to R to improve the condition
% opt.alpha_alg = 'kalman';     % compute alpha using Kalman filtering/smoothing

% 1) run state space GP: O(n)
tic
[post_ss,nlZ_ss,dnlZ_ss] = infState(hyp,mean,cov,lik,x,y,opt);       % inference
[fmu_ss,fs2_ss,ymu_ss,ys2_ss] = post_ss.predict(xs);                % prediction
t_ss = toc;
fprintf('State space GP took %1.2fs on n=%d samples.\n',t_ss,n)

% equivalent computation using the gp function without fast prediction
inf_ss = @(varargin) feval(inf,varargin{:},opt);               % get opt into gp
[nlZ_ss2,dnlZ_ss2,post_ss2] = gp(hyp, inf_ss, mean, {@apxState,cov}, lik, x, y);

if n>2500, fprintf('Large scale; no exact computation possible.\n'), return, end

% 2) run plain dense GP:       O(n^3)
tic
[nlZ,dnlZ,post]   = gp(hyp,inf,mean,cov,lik,x,y);                    % inference
[ymu,ys2,fmu,fs2] = gp(hyp,inf,mean,cov,lik,x,post,xs);             % prediction
t = toc;
fprintf('Plain dense GP took %1.2fs on n=%d samples.\n',t,n)

nor = @(x) norm(x(:)); dev = @(x,z) nor(x-z)/max([1,nor(x),nor(z)]);
vec = @(x) [x.mean; x.lik; x.cov];
err = [dev(post.alpha,post_ss.alpha), dev(post.sW,post_ss.sW),...
       dev(nlZ,nlZ_ss), dev(vec(dnlZ),vec(dnlZ_ss))];

fprintf('err=%1.2e, %1.2e, %1.2e, %1.2e\n',err)

subplot(121), plot(x,y,'b+'), hold on, plot(xs,fmu_ss,'c'), plot(xs,fmu,'b:')
plot(xs,ymu_ss,'m'), plot(xs,ymu,'r:')
legend('training data','f_* state','f_* dense',...
       sprintf('y_* state, logZ=%1.2f',-nlZ),...
       sprintf('y_* dense, logZ=%1.2f',-nlZ_ss))
subplot(122), semilogy(xs,abs(ymu_ss-ymu),'r'), hold on
              semilogy(xs,abs(ys2_ss-ys2),'b')