disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
disp('Hit any key to continue...'); pause

clear all, close all
write_fig = 0;

xg = {linspace(-2,2,120)',linspace(-3,3,120)'};               % construct a grid
[xe,ng] = covGrid('expand',xg);
y = sin(xe(:,2)) + xe(:,1); x = (1:prod(ng))';          % generate training data
xgs = {linspace(-4,4,100)',linspace(-6,6,110)'}; 
y = y + 0.01*gpml_randn(1,size(y));                            % add some noise
[xs,ns] = covGrid('expand',xgs);

cov = {{@covSEiso},{@covSEiso}}; covg = {@covGrid,cov,xg};         % set up a GP
hyp.cov = zeros(4,1); hyp.mean = []; hyp.lik = log(0.1);
hyp = minimize(hyp,@gp,-50,@infGrid,[],covg,[],x,y);  % optimise hyperparameters
opt.cg_maxit = 200; opt.cg_tol = 5e-3;        % perform inference and prediction
post = infGrid(hyp,{@meanZero},covg,'likGauss',x,y,opt); post.L = @(a) a;
ym = gp(hyp, @infGrid,[],covg,[],x,post,xs);

cl = [min(ym(:)), max(ym(:))];
subplot(121), imagesc(xg{1}, xg{2}, reshape(y, ng'),cl), title('data y')
grid on, set(gca,'xtick',-2:2), set(gca,'ytick',-3:3)
subplot(122), imagesc(xgs{1},xgs{2},reshape(ym,ns'),cl), title('prediction ym')
grid on, set(gca,'xtick',-4:4), set(gca,'ytick',-6:6)
rectangle('Position',[-2,-3,4,6])
if write_fig, print -depsc f9.eps; end