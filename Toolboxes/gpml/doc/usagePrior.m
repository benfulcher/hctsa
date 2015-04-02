% demonstrate usage of prior distributions
%
% See also priorDistributions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
clear all, close all

% 1) specify some priors
% a) univariate priors
mu = 1.0; s2 = 0.01^2; nu = 3;
pg = {@priorGauss,mu,s2};                          % Gaussian prior
pl = {'priorLaplace',mu,s2};                        % Laplace prior
pt = {@priorT,mu,s2,nu};                        % Student's t prior
p1 = {@priorSmoothBox1,0,3,15};  % smooth box constraints lin decay
p2 = {@priorSmoothBox2,0,2,15};  % smooth box constraints qua decay
pd = {'priorDelta'}; % fix value of prior exclude from optimisation
pc = {@priorClamped};                         % equivalent to above
lam = 1.05; k = 2.5;
pw = {@priorWeibull,lam,k};                         % Weibull prior

% b) meta priors
pmx = {@priorMix,[0.5,0.5],{pg,pl}};        % mixture of two priors
g = @exp; dg = @exp; ig = @log;
ptr = {@priorTransform,g,dg,ig,pg};    % Gaussian in the exp domain

% c) multivariate priors
m = [1;2]; V = [2,1;1,2];
pG = {@priorGaussMulti,m,V};                    % 2d Gaussian prior
pD = {'priorDeltaMulti'};   % fix value of prior exclude from optim
pC = {@priorClampedMulti};                    % equivalent to above

% 2) evaluation
% pri = pt;   hp = randn(1,3);
% pri = pmx;  hp = randn(1,3);
% pri = ptr;  hp = randn(1,3);
pri = pG;   hp = randn(2,3);

% a) draw a sample from the prior
feval(pri{:})

% b) evaluate prior and derivative if requires
[lp,dlp] = feval(pri{:},hp)

% 3) comprehensive example
x = (0:0.1:10)'; y = 2*x+randn(size(x));   % generate training data
mean = {@meanSum,{@meanConst,@meanLinear}}; % specify mean function
cov = {@covSEiso}; lik = {@likGauss};  % specify covariance and lik
hyp.cov = [log(1);log(1.2)]; hyp.lik = log(0.9); hyp.mean = [2;3];
par = {mean,cov,lik,x,y}; mfun = @minimize; % input for GP function

% a) plain marginal likelihood optimisation (maximum likelihood)
im = @infExact;                                  % inference method
hyp_plain = feval(mfun, hyp, @gp, -10, im, par{:});      % optimise

% b) regularised optimisation (maximum a posteriori) with 1d priors
prior.mean = {pg;pc};  % Gaussian prior for first, clamp second par
prior.cov  = {p1;[]}; % box prior for first, nothing for second par
im = {@infPrior,@infExact,prior};                % inference method
hyp_p1 = feval(mfun, hyp, @gp, -10, im, par{:});         % optimise

% c) regularised optimisation (maximum a posteriori) with Nd priors
prior = [];                                   % clear the structure
% multivariate Student's t prior on the first and second mean hyper
prior.multi{1} = {@priorTMulti,[mu;mu],diag([s2,s2]),nu,...
                 struct('mean',[1,2])};          % use hyper struct
% Equivalent shortcut (same mu and s2 for all dimensions)
prior.multi{1} = {@priorTMulti,mu,s2,nu,struct('mean',[1,2])};
% multivariate Gaussian prior jointly on 1st and 3rd hyper
prior.multi{2} = {@priorGaussMulti,[mu;mu],diag([s2,s2]),...
                 [1,3]};               % use unwrapped hyper vector
% Equivalent shortcut (same mu and s2 for all dimensions)
prior.multi{2} = {@priorGaussMulti,mu,s2,[1,3]};
im = {@infPrior,@infExact,prior};                % inference method
hyp_pN = feval(mfun, hyp, @gp, -10, im, par{:});         % optimise

[unwrap(hyp), unwrap(hyp_plain), unwrap(hyp_p1), unwrap(hyp_pN)]
