% demonstrate usage of classification
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch 2013-10-16.
clear all, close all

%% SAY WHICH CODE WE WISH TO EXERCISE
id = [1,1]; % use Gauss/Exact
id = [1,2; 2,2; 3,2]; % compare Laplace
id = [1,3; 2,3; 3,3]; % study EP
id = [1,5; 2,5]; % look into KL (takes quite a while)
id = [1,4; 2,4; 3,4]; % deal with VB

seed = 943; randn('seed',seed), rand('seed',seed)

ntr = 50; nte = 1e4;                        % number of training and test points
xtr = 10*sort(rand(ntr,1));                                     % sample dataset
p = @(x) 1./(1+exp(-5*sin(x)));                  % "true" underlying probability
ytr = 2*(p(xtr)>rand(ntr,1))-1;                              % draw labels +1/-1
i = randperm(ntr); nout = 3;                                      % add outliers
ytr(i(1:nout)) = -ytr(i(1:nout)); 
xte = linspace(0,10,1e4)';                    % support, we test our function on
cov = {@covSEiso}; sf = 1; ell = 0.7;                             % setup the GP
hyp0.cov  = log([ell;sf]);
mean = {@meanZero};                                                   % m(x) = 0
hyp0.mean = [];
lik_list = {'likGauss','likErf','likLogistic'};          % allowable likelihoods
inf_list = {'infExact','infLaplace','infEP','infVB','infKL'};   % inference algs

Ncg = 50;                                   % number of conjugate gradient steps
sdscale = 0.5;                  % how many sd wide should the error bars become?
col = {'k',[.8,0,0],[0,.5,0],'b',[0,.75,.75],[.7,0,.5]};                % colors
ymu{1} = 2*p(xte)-1; ys2{1} = 0;
for i=1:size(id,1)
  lik = lik_list(id(i,1));                                % setup the likelihood
  if strcmp(lik,'likGauss')
    sn = .2; hyp0.lik = log(sn);
  else
    hyp0.lik = [];
  end
  inf = inf_list{id(i,2)};
  fprintf('OPT: %s/%s\n',lik_list{id(i,1)},inf_list{id(i,2)})
  hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, xtr, ytr);   % opt hypers
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
