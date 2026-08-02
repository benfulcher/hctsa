disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
clear all, close all, write_fig = 0; N = 30;
strc = input(['Which structure?\n  (k)ronecker, (t)oeplitz, (b)ttb, ',...
             '(p)rojection: '],'s');
if isequal(strc,'p'), data = 'r';
else data = input('Which data?\n  (g)rid, (i)ndex, (r)andom: ','s'); end
pred = input('Which prediction?\n  (f)ast, (v)ar+fast, (p)lain: ','s');

% set up the GP
x1 = linspace(-2,2,137)'; x2 = linspace(-3,3,132)';  % construct covariance grid
cov = {{@covSEiso},{@covSEiso}};                % stationary covariance function
mean = {@meanConst}; lik = {@likGauss};     % constant mean, Gaussian likelihood
sf = 1; ell = 0.2; hypcov = log([ell;sf]); hyp.cov = log([ell;sf/2; ell;sf/2]);
hyp.mean = 0.1; sn = 0.1; hyp.lik = log(sn); % mean & likelihood hyperpapameters
switch strc
  case 'k',  xg = {  x1, x2 };                       % plain Kronecker structure
  case 't',  xg = { {x1},x2 };     % Kronecker structure but one Toeplitz factor
  case 'p',  xg = {x2}; cov = {@covSEiso}; hyp.cov = log([ell;sf]); % projection
             hyp.proj = [0,1]; opt.proj = 'ortho';
  otherwise, xg = {{x1,x2}}; cov = {@covSEiso}; hyp.cov = log([ell;sf]);  % BTTB
end
covg = {@apxGrid,cov,xg};                                      % grid covariance
opt.cg_maxit = 500; opt.cg_tol = 1e-5;                          % LCG parameters
inf = @(varargin) infGrid(varargin{:},opt);      % shortcut for inference method

% set up the data
f = @(x) sin(x(:,2)) + x(:,1);                        % 2d function to be learnt
xe = apxGrid('expand',xg); idx = randperm(size(xe,1))'; idx = idx(1:fix(end/5));
switch data
  case 'g'                       % located on covariance grid but 2d coordinates
    x = xe(idx,:);                     xx = x;         sdata = 'gridded';
  case 'i'                              % located on covariance grid but indices
    x = idx;                           xx = xe(idx,:); sdata = 'grid index';
  otherwise                                             % arbitrary 2d locations
    x = (2*rand(2e4,2)-1)*diag([2,3]); xx = x;         sdata = 'random';
end
fprintf('Use %s training data.\n',sdata)
y = f(xx); y = y + 0.1*randn(size(y));  % add some observation noise to the data

% construct a grid covering the training data from scratch
if isequal(strc,'k'), xg = apxGrid('create',xx,true,[137,132]); end

% set up the query
xgs = {linspace(-4,4,400)',linspace(-6,6,410)'};    % construct a test data grid
[xs,ns] = apxGrid('expand',xgs);
par = {mean,covg,lik,x};              % shortcut for Gaussian process parameters
opt.stat = true;                   % show some more information during inference
opt.ndcovs = 25;                    % ask for sampling-based (exact) derivatives

% if strcmp(pred,'v'), opt.pred_var = -50; end % Lanczos-based estimation (LOVE)
if strcmp(pred,'v'), opt.pred_var =  20; end % sampling-based estimation

tic, [post,nlZ,dnlZ] = infGrid(hyp,par{:},y,opt); ti = toc; tic  % run inference
switch pred
  case 'f',                     [fmu,fs2,ymu,ys2] = post.predict(xs);
  case 'v',                     [fmu,fs2,ymu,ys2] = post.predict(xs);
  otherwise, post.L = @(a) 0*a; [ymu,ys2,fmu,fs2] = gp(hyp,inf,par{:},post,xs);
end
fprintf('Inference/prediction took %1.2f/%1.2f[s]\n',ti,toc)

clm = [min(ymu(:)), max(ymu(:))];
cls = [sn, sqrt(sf^2+sn^2)];
scrsz = get(0,'ScreenSize');
figure('Position',[200 200 scrsz(3)/2 scrsz(4)/2])
subplot(131), scatter(xx(:,1),xx(:,2),3,y,'filled','s')
 title(sprintf('%s training data y, n=%d',sdata,numel(y)))
 set(gca,'xtick',-2:2), set(gca,'ytick',-3:3)
 rectangle('Position',[-2,-3,4,6],'EdgeColor','r','LineWidth',3)
subplot(132), imagesc(xgs{1},xgs{2},reshape(ymu,ns')',clm)
 title(sprintf('test set prediction mean, n_*=%d',numel(ymu)))
 axis xy, grid on, set(gca,'xtick',-4:4), set(gca,'ytick',-6:6)
 rectangle('Position',[-2,-3,4,6],'EdgeColor','r','LineWidth',3)
subplot(133), imagesc(xgs{1},xgs{2},reshape(sqrt(ys2),ns')',cls)
 title(sprintf('test set prediction std, n_*=%d',numel(ymu)))
 axis xy, grid on, set(gca,'xtick',-4:4), set(gca,'ytick',-6:6)
 rectangle('Position',[-2,-3,4,6],'EdgeColor','r','LineWidth',3)
if write_fig, print -depsc f9.eps; end