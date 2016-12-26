function out = MF_hmm_fit(y,trainp,numStates,randomSeed)
% MF_hmm_fit    Fits a Hidden Markov Model to sequential data.
%
% Actually highly stochastic, so for reproducible results helps to reset the
% random seed...
%
%---INPUTS:
% y, the input time series
% trainp, the proportion of data to train on, 0 < trainp < 1
% numStates, the number of states in the HMM

% Uses Zoubin Gharamani's implementation of HMMs for real-valued Gaussian
% observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz
%
% Uses ZG_hmm (renamed from hmm) and ZG_hmm_cl (renamed from hmm_cl)
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

% Check required function files exist:
if ~exist(fullfile('ZG_hmm','ZG_hmm_cl'),'file') || ~exist(fullfile('ZG_hmm','ZG_hmm_cl'),'file')
    error('Could not find the required HMM fitting functions (Zoubin Gharamani''s code)');
end

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % number of samples in time series

if nargin < 2 || isempty(trainp)
    fprintf(1,'Training on 80%% of the data by default\n');
    trainp = 0.8; % train on 80% of the data
end

if nargin < 3 || isempty(numStates)
    fprintf(1,'Using 3 states by default\n');
    numStates = 3; % use 3 states
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
Ntrain = floor(trainp*N);
if Ntrain == N
    error('No data for test set for HMM fitting?!')
end

ytrain = y(1:Ntrain);
if Ntrain < N
    ytest = y(Ntrain+1:end);
    Ntest = length(ytest);
end

% Train HMM with <numStates> states for 30 cycles of EM (or until
% convergence); default termination tolerance
[Mu, Cov, P, Pi, LL] = ZG_hmm(ytrain,Ntrain,numStates,30);

% ------------------------------------------------------------------------------
%% Output statistics on the training
% ------------------------------------------------------------------------------

% Mean vector, Mu
Musort = sort(Mu,'ascend');
for i = 1:length(Mu)
    out.(sprintf('Mu_%u',i)) = Musort(i); % use dynamic field referencing
    % eval(sprintf('out.Mu_%u = Musort(%u);',i,i));
end
out.meanMu = mean(Mu);
out.rangeMu = max(Mu) - min(Mu);
out.maxMu = max(Mu);
out.minMu = min(Mu);

% Covariance Cov
out.Cov = Cov;

% Transition matrix
out.Pmeandiag = mean(diag(P));
out.stdmeanP = std(mean(P));
out.maxP = max(P(:));
out.meanP = mean(P(:)); % I guess this is just 1/numStates? A constant so not a useful output.
out.stdP = std(P(:));

% Within-sample log-likelihood
out.LLtrainpersample = max(LL)/Ntrain; % loglikelihood per sample
out.nit = length(LL); % number of iterations

% ------------------------------------------------------------------------------
%% Calculate log likelihood for the test data
% ------------------------------------------------------------------------------
lik = ZG_hmm_cl(ytest,Ntest,numStates,Mu,Cov,P,Pi);

out.LLtestpersample = lik/Ntest;

out.LLdifference = out.LLtestpersample - out.LLtrainpersample;


end
