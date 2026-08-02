% demonstrate usage of regression
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2016-08-29.
clear all, close all

%% SAY WHICH CODE WE WISH TO EXERCISE
id = [1,1]; % use Gauss/Exact
id = [1,2; 3,2; 4,2]; % compare Laplace
id = [1,3; 2,3; 3,3]; % study EP
id = [1,5; 2,5]; % look into KL (takes quite a while)
id = [1,4; 2,4; 3,4; 4,4]; % deal with VB

seed = 197; randn('seed',seed), rand('seed',seed)

ntr = 50; nte = 1e4;                        % number of training and test points
xtr = 10*sort(rand(ntr,1));                                     % sample dataset
f = @(x) sin(x)+sqrt(x);                            % "true" underlying function
sn = 0.2;
ytr = f(xtr) + randn(ntr,1)*sn;                             % add Gaussian noise
i = randperm(ntr); nout = 3;                                      % add outliers
ytr(i(1:nout)) = 5;
xte = linspace(0,10,1e4)';                    % support, we test our function on

cov = {@covSEiso}; sf = 1; ell = 0.4;                             % setup the GP
hyp0.cov  = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; a = 1/5; b = 1;       % m(x) = a*x+b
hyp0.mean = [a;b];

lik_list = {'likGauss','likLaplace','likSech2','likT'};   % possible likelihoods
inf_list = {'infGaussLik','infLaplace','infEP','infVB','infKL'};% inference algs

Ncg = 50;                                   % number of conjugate gradient steps
sdscale = 0.5;                  % how many sd wide should the error bars become?
col = {'k',[.8,0,0],[0,.5,0],'b',[0,.75,.75],[.7,0,.5]};                % colors
ymu{1} = f(xte); ys2{1} = sn^2; nlZ(1) = -Inf;
for i=1:size(id,1)
  lik = lik_list{id(i,1)};                                % setup the likelihood
  if strcmp(lik,'likT')
    nu = 4;
    hyp0.lik  = log([nu-1;sqrt((nu-2)/nu)*sn]);
  else
    hyp0.lik  = log(sn);
  end
  inf = inf_list{id(i,2)};
  fprintf('OPT: %s/%s\n',lik_list{id(i,1)},inf_list{id(i,2)})
  if Ncg==0
    hyp = hyp0;
  else
    hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtr, ytr); % opt hypers
  end
  [ymu{i+1}, ys2{i+1}] = gp(hyp, inf, mean, cov, lik, xtr, ytr, xte);  % predict
  [nlZ(i+1)] = gp(hyp, inf, mean, cov, lik, xtr, ytr);
end

figure, hold on
for i=1:size(id,1)+1
  plot(xte,ymu{i},'Color',col{i},'LineWidth',2)
  if i==1
    leg = {'function'};
  else
    leg{end+1} = sprintf('%s/%s -lZ=%1.2f',...
                                lik_list{id(i-1,1)},inf_list{id(i-1,2)},nlZ(i));
  end
end
for i=1:size(id,1)+1
  ysd = sdscale*sqrt(ys2{i});
  fill([xte;flipud(xte)],[ymu{i}+ysd;flipud(ymu{i}-ysd)],...
       col{i},'EdgeColor',col{i},'FaceAlpha',0.1,'EdgeAlpha',0.3);
end
for i=1:size(id,1)+1, plot(xte,ymu{i},'Color',col{i},'LineWidth',2), end
plot(xtr,ytr,'k+'), plot(xtr,ytr,'ko'), legend(leg)