disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
disp('Hit any key to continue...'); pause

disp(' '); disp('clear all, close all')
clear all, close all
write_fig = 0;
disp(' ')

disp('n1 = 80; n2 = 40;                   % number of data points from each class')
n1 = 80; n2 = 40;
disp('S1 = eye(2); S2 = [1 0.95; 0.95 1];           % the two covariance matrices')
S1 = eye(2); S2 = [1 0.95; 0.95 1];
disp('m1 = [0.75; 0]; m2 = [-0.75; 0];                            % the two means')
m1 = [0.75; 0]; m2 = [-0.75; 0];
disp(' ')

disp('x1 = bsxfun(@plus, chol(S1)''*gpml_randn(0.2, 2, n1), m1);')
x1 = bsxfun(@plus, chol(S1)'*gpml_randn(0.2, 2, n1), m1);
disp('x2 = bsxfun(@plus, chol(S2)''*gpml_randn(0.3, 2, n2), m2);')         
x2 = bsxfun(@plus, chol(S2)'*gpml_randn(0.3, 2, n2), m2);         
disp(' ')

disp('x = [x1 x2]''; y = [-ones(1,n1) ones(1,n2)]'';')
x = [x1 x2]'; y = [-ones(1,n1) ones(1,n2)]';
figure(6)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on;');
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'');');
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12);
disp(' ')

disp('[t1 t2] = meshgrid(-4:0.1:4,-4:0.1:4);')
[t1 t2] = meshgrid(-4:0.1:4,-4:0.1:4);
disp('t = [t1(:) t2(:)]; n = length(t);               % these are the test inputs')
t = [t1(:) t2(:)]; n = length(t);               % these are the test inputs
disp('tmm = bsxfun(@minus, t, m1'');')
tmm = bsxfun(@minus, t, m1');
disp('p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));')
p1 = n1*exp(-sum(tmm*inv(S1).*tmm/2,2))/sqrt(det(S1));
disp('tmm = bsxfun(@minus, t, m2'');')
tmm = bsxfun(@minus, t, m2');
disp('p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));')
p2 = n2*exp(-sum(tmm*inv(S2).*tmm/2,2))/sqrt(det(S2));
set(gca, 'FontSize', 24)
disp('contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])')
contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.1:0.1:0.9])
[c h] = contour(t1, t2, reshape(p2./(p1+p2), size(t1)), [0.5 0.5]);
set(h, 'LineWidth', 2)
colorbar
grid
axis([-4 4 -4 4])
if write_fig, print -depsc f6.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp(' ')
disp('meanfunc = @meanConst; hyp.mean = 0;')
meanfunc = @meanConst; hyp.mean = 0;
disp('covfunc = @covSEard;   hyp.cov = log([1 1 1]);')
covfunc = @covSEard;   hyp.cov = log([1 1 1]);
disp('likfunc = @likErf;')
likfunc = @likErf;
disp(' ')

disp('hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);')
hyp = minimize(hyp, @gp, -40, @infEP, meanfunc, covfunc, likfunc, x, y);
disp('[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n, 1));')
[a b c d lp] = gp(hyp, @infEP, meanfunc, covfunc, likfunc, x, y, t, ones(n,1));
disp(' ')
figure(7)
set(gca, 'FontSize', 24)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on')
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'')')
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
[c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
set(h, 'LineWidth', 2)
colorbar
grid
axis([-4 4 -4 4])
if write_fig, print -depsc f7.eps; end
disp(' '); disp('Hit any key to continue...'); pause

disp('large scale classification using the FITC approximation')
disp('[u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)];')
[u1,u2] = meshgrid(linspace(-2,2,5)); u = [u1(:),u2(:)]; clear u1; clear u2
disp('nu = size(u,1);')
nu = size(u,1);
disp('covfuncF = {@covFITC, {covfunc}, u};')
covfuncF = {@covFITC, {covfunc}, u};
disp('inffunc = @infFITC_Laplace;')
inffunc = @infFITC_EP;                     % one could also use @infFITC_Laplace
disp('hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfuncF, likfunc, x, y);')
hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfuncF, likfunc, x, y);
disp('[a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));')
[a b c d lp] = gp(hyp, inffunc, meanfunc, covfuncF, likfunc, x, y, t, ones(n,1));
disp(' ')
figure(8)
set(gca, 'FontSize', 24)
disp('plot(x1(1,:), x1(2,:), ''b+''); hold on')
plot(x1(1,:), x1(2,:), 'b+', 'MarkerSize', 12); hold on
disp('plot(x2(1,:), x2(2,:), ''r+'')')
plot(x2(1,:), x2(2,:), 'r+', 'MarkerSize', 12)
disp('contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)')
contour(t1, t2, reshape(exp(lp), size(t1)), 0.1:0.1:0.9)
[c h] = contour(t1, t2, reshape(exp(lp), size(t1)), [0.5 0.5]);
set(h, 'LineWidth', 2)
disp('plot(u(:,1),u(:,2),''ko'', ''MarkerSize'', 12)')
plot(u(:,1),u(:,2),'ko', 'MarkerSize', 12)
colorbar
grid
axis([-4 4 -4 4])
if write_fig, print -depsc f8.eps; end
disp(' ')