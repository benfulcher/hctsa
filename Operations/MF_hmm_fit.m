function out = MF_hmm(y,trainp,nstates)
% Uses Zoubin Gharamani's implementation of HMMs for real-valued Gaussian
% observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz

% I've wrapped it up to work and give a structure as output for
% implementation in the time series database...
% Ben Fulcher 9/4/2010


%% Check Inputs
N = length(y); % number of samples in time series

if nargin < 2 || isempty(trainp)
    fprintf(1,'Training on 80%% of the data by default\n')
    trainp = 0.8; % train on 80% of the data
end

if nargin < 3 || isempty(nstates)
    fprintf(1,'Using 3 states by default\n')
    nstates = 3; % use 3 states
end

%% Train the HMM
% divide up dataset into training (ytrain) and test (ytest) portions
Ntrain = floor(trainp*N);
if Ntrain == N
    error('No data for test set for HMM fitting?!')
end

ytrain = y(1:Ntrain);
if Ntrain < N
    ytest = y(Ntrain+1:end);
    Ntest = length(ytest);
end


% train HMM with <nstates> states for 30 cycles of EM (or until
% convergence); default termination tolerance
[Mu, Cov, P, Pi, LL] = ZG_hmm(ytrain,Ntrain,nstates,30);

%% Output statistics on the training

% mean vector, Mu
Musort = sort(Mu,'ascend');
for i = 1:length(Mu)
    eval(sprintf('out.Mu_%u = Musort(%u);',i,i));
end
out.meanMu = mean(Mu);
out.rangeMu = max(Mu) - min(Mu);
% out.maxMu = max(Mu);
% out.minMu = min(Mu);


% Covariance Cov
out.Cov = Cov;

% transition matrix
out.Pmeandiag = mean(diag(P));
out.stdmeanP = std(mean(P));
out.maxP = max(P(:));
out.meanP = mean(P(:));
out.stdP = std(P(:));

% Within-sample log-likelihood
out.LLtrainpersample = max(LL)/Ntrain; % loglikelihood per sample
out.nit = length(LL); % number of iterations

%% Calculate log likelihood for the test data
lik = ZG_hmm_cl(ytest,Ntest,nstates,Mu,Cov,P,Pi);

out.LLtestpersample = lik/Ntest;

out.LLdifference = out.LLtestpersample - out.LLtrainpersample;


end