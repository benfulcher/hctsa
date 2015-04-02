% demonstrate usage of likelihood functions
%
% See also likFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
clear all, close all
n = 5; f = randn(n,1);       % create random latent function values

% set up simple classification likelihood functions
yc = sign(f);
lc0 = {'likErf'};     hypc0 = [];   % no hyperparameters are needed
lc1 = {@likLogistic}; hypc1 = [];    % also function handles are OK
lc2 = {'likUni'};     hypc2 = [];
lc3 = {'likMix',{'likUni',@likErf}}; hypc3 = log([1;2]); %mixture

% set up simple regression likelihood functions
yr = f + randn(n,1)/20;
sn = 0.1;                                % noise standard deviation
lr0 = {'likGauss'};   hypr0 = log(sn);
lr1 = {'likLaplace'}; hypr1 = log(sn);
lr2 = {'likSech2'};   hypr2 = log(sn);
nu = 4;                              % number of degrees of freedom
lr3 = {'likT'};       hypr3 = [log(nu-1); log(sn)];
lr4 = {'likMix',{lr0,lr1}}; hypr4 = [log([1,2]);hypr0;hypr1];

a = 1; % set up warped Gaussian with g(y) = y + a*sign(y).*y.^2
lr5 = {'likGaussWarp',['poly2']}; hypr5 = log([a;sn]);
lr6 = {'likGumbel','+'}; hypr6 = log(sn);

% set up Poisson regression
yp = fix(abs(f)) + 1;
lp0 = {@likPoisson,'logistic'}; hypp0 = [];
lp1 = {@likPoisson,'exp'};      hypp1 = [];

% set up other GLM likelihoods for positive or interval regression
lg1 = {@likGamma,'logistic'}; al = 2;    hyp.lik = log(al);
lg2 = {@likInvGauss,'exp'};   lam = 1.1; hyp.lik = log(lam);
lg3 = {@likBeta,'expexp'};    phi = 2.1; hyp.lik = log(phi);
lg4 = {@likBeta,'logit'};     phi = 4.7; hyp.lik = log(phi);

% 0) specify the likelihood function
lik = lc0; hyp = hypc0; y = yc;
% lik = lr4; hyp = hypr4; y = yr;
% lik = lp1; hyp = hypp1; y = yp;

% 1) query the number of parameters
feval(lik{:})

% 2) evaluate the likelihood function on f
exp(feval(lik{:},hyp,y,f))

% 3a) evaluate derivatives of the likelihood
[lp,dlp,d2lp,d3lp] = feval(lik{:}, hyp, y, f, [], 'infLaplace');

% 3b) compute Gaussian integrals w.r.t. likelihood
mu = f; s2 = rand(n,1);
[lZ,dlZ,d2lZ] = feval(lik{:}, hyp, y, mu, s2, 'infEP');

% 3c) obtain lower bound on likelihood
ga = rand(n,1);
[b,z] = feval(lik{:}, hyp, y, [], ga, 'infVB');
