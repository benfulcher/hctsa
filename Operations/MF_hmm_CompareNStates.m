function out = MF_hmm_CompareNStates(y,trainp,nstater,randomSeed)
% MF_hmm_CompareNStates     Hidden Markov Model (HMM) fitting to a time series.
%
% Fits HMMs with different numbers of states, and compares the resulting
% test-set likelihoods.
%
% The code relies on Zoubin Gharamani's implementation of HMMs for real-valued
% Gassian-distributed observations, including the hmm and hmm_cl routines (
% renamed ZG_hmm and ZG_hmm_cl here).
% Implementation of HMMs for real-valued Gaussian observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz
%
%---INPUTS:
%
% y, the input time series
%
% trainp, the initial proportion of the time series to train the model on
%
% nstater, the vector of state numbers to compare. E.g., (2:4) compares a number
%               of states 2, 3, and 4.
%
%---OUTPUTS: statistics on how the log likelihood of the test data changes with
% the number of states n_{states}$. We implement the code for p_{train} = 0.6$
% as n_{states}$ varies across the range n_{states} = 2, 3, 4$.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % number of samples in time series

if nargin < 2 || isempty(trainp)
    fprintf(1,'Training the model on 60%% of the data by default\n');
    trainp = 0.6; % train on 60% of the data
end
Ntrain = floor(trainp*N); % number of initial samples to train the model on

if nargin < 3 || isempty(nstater)
    fprintf(1,'Using 2--4 states by default\n');
    nstater = (2:4); % use 2:4 states
end

if nargin < 4
    randomSeed = [];
end

%-------------------------------------------------------------------------------
% Deal with random seeds
%-------------------------------------------------------------------------------
BF_ResetSeed(randomSeed); % reset the random seed if specified

% ------------------------------------------------------------------------------
%% Train the HMM
% ------------------------------------------------------------------------------
% Divide up dataset into training (ytrain) and test (ytest) portions
ytrain = y(1:Ntrain);
if Ntrain < N
    ytest = y(Ntrain+1:end);
    Ntest = length(ytest);
end

Nstate = length(nstater);
LLtrains = zeros(Nstate,1);
LLtests = zeros(Nstate,1);

for j = 1:Nstate
    numStates = nstater(j);
    % train HMM with <numStates> states for 30 cycles of EM (or until
    % convergence); default termination tolerance
    [Mu, Cov, P, Pi, LL] = ZG_hmm(ytrain,Ntrain,numStates,30);

    LLtrains(j) = LL(end)/Ntrain;

    %% Calculate log likelihood for the test data
    lik = ZG_hmm_cl(ytest,Ntest,numStates,Mu,Cov,P,Pi);

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

for i = 1:Nstate-1
    out.(sprintf('LLtestdiff%u',i)) = LLtests(i+1) - LLtests(i);
end


end
