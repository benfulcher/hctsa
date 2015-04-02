% demonstrate usage of mean functions
%
% See also meanFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
clear all, close all
n = 5; D = 2; x = randn(n,D);            % create a random data set

% set up simple mean functions
m0 = {'meanZero'};  hyp0 = [];      % no hyperparameters are needed
m1 = {'meanOne'};   hyp1 = [];      % no hyperparameters are needed
mc = {@meanConst};  hypc = 2;  % also function handles are possible
ml = {@meanLinear}; hypl = [2;3];              % m(x) = 2*x1 + 3*x2
mp = {@meanPoly,2}; hypp = [1;1;2;3];  % m(x) = x1+x2+2*x1^2+3*x2^2
mn = {@meanNN,[1,0; 0,1],[0.9,0.5]}; hypn = [];  % nearest neighbor
s = 12; hypd = randn(s,1);           % discrete mean with 12 hypers
md = {'meanDiscrete',s};
hyp.cov = [0;0]; hypg = [];                    % GP predictive mean
xt = randn(2*n,D); yt = sign(xt(:,1)-xt(:,2));      % training data
mg = {@meanGP,hyp,@infEP,@meanZero,@covSEiso,@likErf,xt,yt};
hype = [0;0; log(0.1)];             % regression GP predictive mean
xt = randn(2*n,D); yt = xt(:,1).*xt(:,2);           % training data
me = {@meanGPexact,@meanZero,@covSEiso,xt,yt};

% set up composite mean functions
msc = {'meanScale',{m1}};      hypsc = [3; hyp1];      % scale by 3
msu = {'meanSum',{m0,mc,ml}};  hypsu = [hyp0; hypc; hypl];    % sum
mpr = {@meanProd,{mc,ml}};     hyppr = [hypc; hypl];      % product
mpo = {'meanPow',3,msu};       hyppo = hypsu;         % third power
mask = [false,true];     % mask excluding all but the 2nd component
mma = {'meanMask',mask,ml};    hypma = hypl(mask);
mpf = {@meanPref,ml};          hyppf = 2;  % linear pref with slope

% 0) specify mean function
% mean = md; hyp = hypd; x = randi([1,s],n,1);
% mean = mn; hyp = hypn;
% mean = mg; hyp = hypg;
mean = me; hyp = hype;
% mean = m0;  hyp = hyp0;
% mean = msu; hyp = hypsu;
% mean = mpr; hyp = hyppr;
% mean = mpo; hyp = hyppo;
% mean = mpf; hyp = hyppf;

% 1) query the number of parameters
feval(mean{:})

% 2) evaluate the function on x
feval(mean{:},hyp,x)

% 3) compute the derivatives w.r.t. to hyperparameter i
i = 2; feval(mean{:},hyp,x,i)
