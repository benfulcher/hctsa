clear all, close all
n1=10; n2=5;                       % number of data points from each class n1/n2
sc=1;                                  % change the scale of the initial problem

% 1) draw data from two Gaussians
S1 = eye(2); S2 = [1 0.95; 0.95 1];
m1 = [0.75;  0]; m2 = [-0.75; 0];
randn('seed',17); rand('seed',17)
% coordinates
x1 = chol(S1)'*randn(2,n1)+repmat(m1,1,n1);
x2 = chol(S2)'*randn(2,n2)+repmat(m2,1,n2);
x1(:,5) = x1(:,2); x1(:,1) = x1(:,4); x2(:,1) = x2(:,2);  % rank deficiency in K
xtr = [x1 x2]';
% class labels
ytr = [repmat(-1,[1,n1]) repmat(1,[1,n2])]';
tt = (-5:.3:5)*sc; [t1 t2] = meshgrid(tt,tt);
x1 = x1*sc; x2=x2*sc; xtr=xtr*sc; m1=m1*sc; m2=m2*sc;
t  = [t1(:) t2(:)];
z1 = exp(-sum((t-repmat(m1',length(t),1))*inv(S1).*(t-repmat(m1',length(t),1)),2)/2);
z2 = exp(-sum((t-repmat(m2',length(t),1))*inv(S2).*(t-repmat(m2',length(t),1)),2)/2);
xte = t; clear t
xtr(end,:)=[2.9,2.3]; ytr(end)=-1;                                       % modif

% 2) set GP parameters
cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
lik =  {'likLogistic'};  hyp.lik  = [];                  % likLogistic or likErf
mn  = {'meanZero'};      hyp.mean = [];
inf = 'infEP';

% 3) predict using approximated inference
[nlZ,dnlZ,post] = gp(hyp, inf, mn, cov, lik, xtr, ytr);
[ymu,ys2,fmu,fs2,junk,post] = gp(hyp,inf,mn,cov,lik,xtr,post,xte);

% 4) prediction using MCMC sampling
% set MCMC parameters, see some more details in inf/infMCMC.m
% We have two samplers implemented, namely
%  hmc - Hybrid Monte Carlo, and
%  ess - Elliptical Slice Sampling.
par.sampler = 'hmc'; par.Nsample = 20;
% par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
tic
[posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,xtr,ytr,par);
toc
[ymus,ys2s,fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,xtr,posts,xte);

% 5a) echo results
fprintf('nlZ-EP=%f, nlZ-AIS=%f\n', nlZ, nlZs)
fprintf('acceptance rate (MCMC) = %1.2f%%\n',100*posts.acceptance_rate_MCMC)
for r=1:length(posts.acceptance_rate_AIS)
  fprintf('acceptance rate (AIS) = %1.2f%%\n',100*posts.acceptance_rate_AIS(r))
end
% 5b) print results
figure
subplot(221)
  plot([-12,12],[-12,12],'r'), hold on, plot(fmus,fmu,'k.'), title('\mu_f')
  xlabel('MCMC'), ylabel(inf)
subplot(222)
  plot([0,10],[0,10],'r'), hold on, plot(sqrt(fs2s),sqrt(fs2),'k.')
  title('\sigma_f')
  xlabel('MCMC'), ylabel(inf)
subplot(223)
  plot([-1,1],[-1,1],'r'), hold on, plot(ymus,ymu,'k.'), title('\mu_y')
  xlabel('MCMC'), ylabel(inf)
subplot(224)
  plot([0,1],[0,1],'r'), hold on, plot(sqrt(ys2s),sqrt(ys2),'k.')
  title('\sigma_y')
  xlabel('MCMC'), ylabel(inf)