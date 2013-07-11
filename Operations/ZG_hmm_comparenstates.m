function out = ZG_hmm_comparenstates(y,trainp,nstater)
% Uses Zoubin Gharamani's implementation of HMMs for real-valued Gaussian
% observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz

% I've wrapped it up to work and give a structure as output for
% implementation in the time series database...

% This routine fits hmms with a different number of states and compares the
% test-set likelihoods.

% Ben Fulcher 9/4/2010


%% Check Inputs
N = length(y); % number of samples in time series

if nargin < 2 || isempty(trainp)
    fprintf(1,'Training the model on 60% of the data by default\n')
    trainp = 0.6; % train on 60% of the data
end

if nargin < 3 || isempty(nstater)
    fprintf(1,'Using 2--4 states by default\n')
    nstater = (2:4); % use 2:4 states
end


%% Train the HMM
% divide up dataset into training (ytrain) and test (ytest) portions
Ntrain = floor(trainp*N);
ytrain = y(1:Ntrain);
if Ntrain < N
    ytest = y(Ntrain+1:end);
    Ntest = length(ytest);
end

Nstate = length(nstater);
LLtrains = zeros(Nstate,1);
LLtests = zeros(Nstate,1);

for j = 1:Nstate
    nstates = nstater(j);
    % train HMM with <nstates> states for 30 cycles of EM (or until
    % convergence); default termination tolerance
    [Mu, Cov, P, Pi, LL] = hmm(ytrain,Ntrain,nstates,30);
    
    LLtrains(j) = LL(end)/Ntrain;
    
    %% Calculate log likelihood for the test data
    lik = hmm_cl(ytest,Ntest,nstates,Mu,Cov,P,Pi);

    LLtests(j) = lik/Ntest;
end

%% Output some statistics

out.meanLLtrain = mean(LLtrains);
out.meanLLtest = mean(LLtests);
out.maxLLtrain = max(LLtrains);
out.maxLLtest = max(LLtests);
out.chLLtrain = LLtrains(end)-LLtrains(1);
out.chLLtest = LLtests(end)-LLtests(1);
out.meandiffLLtt = mean(abs(LLtests-LLtrains));

for i = 1 : Nstate-1
    eval(sprintf('out.LLtestdiff%u = LLtests(%u) - LLtests(%u);',i,i+1,i));
end


end