% demonstrate usage of mean functions
%
% See also meanFunctions.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2013-01-21
clear all, close all
n = 5; D = 2; x = randn(n,D);            % create a random data set

% set up simple mean functions
m0 = {'meanZero'};  hyp0 = [];      % no hyperparameters are needed
m1 = {'meanOne'};   hyp1 = [];      % no hyperparameters are needed
mc = {@meanConst};  hypc = 2;  % also function handles are possible
ml = {@meanLinear}; hypl = [2;3];              % m(x) = 2*x1 + 3*x2

% set up composite mean functions
msc = {'meanScale',{m1}};      hypsc = [3; hyp1];      % scale by 3
msu = {'meanSum',{m0,mc,ml}};  hypsu = [hyp0; hypc; hypl];    % sum
mpr = {@meanProd,{mc,ml}};     hyppr = [hypc; hypl];      % product
mpo = {'meanPow',{3,msu}};     hyppo = hypsu;         % third power
mask = [0,1,0]; %   binary mask excluding all but the 2nd component
mma = {'meanMask',{mask,mpo{:}}}; hypma = hyppo;

% 0) specify mean function
% mean = m0;  hyp = hyp0;
% mean = msu; hyp = hypsu;
% mean = mpr; hyp = hyppr;
mean = mpo; hyp = hyppo;

% 1) query the number of parameters
feval(mean{:})

% 2) evaluate the function on x
feval(mean{:},hyp,x)

% 3) compute the derivatives w.r.t. to hyperparameter i
i = 2; feval(mean{:},hyp,x,i)
