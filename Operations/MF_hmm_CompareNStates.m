function out = MF_hmm_CompareNStates(y, trainp, nstater)
% MF_hmm_CompareNStates   How the fit of hidden Markov models to the series changes with the number of hidden states.
%
% Fits Gaussian hidden Markov models (HMMs) with different numbers of states to
% the first trainp proportion of the time series (each with at most 30 cycles of
% EM), and compares the resulting log-likelihoods per sample on the training part
% and on the held-out remainder. Each model is fitted deterministically, by the best
% of six fixed starting points with a floor on the shared variance (ZG_hmm_fit, as
% in MF_hmm_Fit).
%
% The code relies on Zoubin Gharamani's implementation of HMMs for real-valued
% Gaussian-distributed observations, including the hmm and hmm_cl routines (
% renamed ZG_hmm and ZG_hmm_cl here), and ZG_hmm_fit.
% Implementation of HMMs for real-valued Gaussian observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz
%
% ---INPUTS:
%
% y, the input time series
%
% trainp, the initial proportion of the time series to train the model on
%         (default: 0.6)
%
% nstater, the vector of numbers of states to compare (default: 2:4)
%
% ---OUTPUTS:
% meanLLtrain, maxLLtrain: mean and maximum across models of the log-likelihood per
%       sample on the training part
% meanLLtest, maxLLtest: mean and maximum across models of the log-likelihood per
%       sample on the test part
% chLLtrain, chLLtest: change in training and test log-likelihood per sample from
%       the model with the fewest states to the one with the most
% meandiffLLtt: mean across models of the absolute difference between the test and
%       training log-likelihoods per sample
% LLtestdiff1, LLtestdiff2, ...: change in test log-likelihood per sample from the
%       i-th to the (i+1)-th number of states in nstater (one fewer than the
%       number of models)

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
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
	fprintf(1, 'Training the model on 60%% of the data by default\n');
	trainp = 0.6; % train on 60% of the data
end
Ntrain = floor(trainp * N); % number of initial samples to train the model on

if nargin < 3 || isempty(nstater)
	fprintf(1, 'Using 2--4 states by default\n');
	nstater = (2:4); % use 2:4 states
end

% ------------------------------------------------------------------------------
%% Train the HMM
% ------------------------------------------------------------------------------
% Divide up dataset into training (yTrain) and test (yTest) portions
if Ntrain >= N
	% Every subsequent step evaluates the fitted models on a held-out test
	% portion, so a training proportion that consumes the whole series leaves
	% nothing to evaluate on (this used to leave yTest/Ntest undefined and
	% error deeper in the loop):
	error('trainp = %g leaves no test data for a series of length %u', trainp, N);
end
if Ntrain < 2
	% Data-dependent: too short to fit anything on the training portion
	warning('Time series (N = %u) too short to train on %g of it', N, trainp);
	out = NaN; return
end
yTrain = y(1:Ntrain);
yTest = y(Ntrain + 1:end);
Ntest = length(yTest);

Nstate = length(nstater);
LLtrains = zeros(Nstate, 1);
LLtests = zeros(Nstate, 1);

for j = 1:Nstate
	numStates = nstater(j);
	% train HMM with <numStates> states for 30 cycles of EM (or until
	% convergence), from fixed starting points; default termination tolerance
	[Mu, Cov, P, Pi, LL] = ZG_hmm_fit(yTrain, numStates, 30);

	LLtrains(j) = LL(end) / Ntrain;

	%% Calculate log likelihood for the test data
	lik = ZG_hmm_cl(yTest, Ntest, numStates, Mu, Cov, P, Pi);

	LLtests(j) = lik / Ntest;
end

%% Output some statistics
out.meanLLtrain = mean(LLtrains);
out.meanLLtest = mean(LLtests);
out.maxLLtrain = max(LLtrains);
out.maxLLtest = max(LLtests);
out.chLLtrain = LLtrains(end) - LLtrains(1);
out.chLLtest = LLtests(end) - LLtests(1);
out.meandiffLLtt = mean(abs(LLtests - LLtrains));

for i = 1:Nstate - 1
	out.(sprintf('LLtestdiff%u', i)) = LLtests(i + 1) - LLtests(i);
end

end
