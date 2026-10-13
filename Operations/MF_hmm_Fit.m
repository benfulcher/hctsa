function out = MF_hmm_Fit(y, trainp, numStates)
% MF_hmm_Fit   A hidden Markov model fitted to the first part of the series, and how well it describes the rest.
%
% Fits a hidden Markov model (HMM) with Gaussian emissions to the first trainp
% proportion of the time series, then evaluates its likelihood on the remainder.
% The emissions of all states share one variance (a tied covariance). The model is
% trained with at most 30 cycles of EM (Baum-Welch), or until convergence.
%
% The fit is deterministic. EM is run from six fixed starting points (state means at
% the quantiles (k-1/2)/numStates of the training data, or evenly spaced within one
% standard deviation of its mean; a variance equal to that of the training data; and
% probabilities 0.5, 0.9 and 0.99 of staying in a state); the fit with the highest
% training log-likelihood is kept (ZG_hmm_fit). The
% shared variance cannot fall below 1% of the variance of the training data, so that
% states placed on a few repeated values do not give unbounded likelihoods.
%
% Uses Zoubin Gharamani's implementation of HMMs for real-valued Gaussian
% observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz
%
% Uses ZG_hmm (renamed from hmm), ZG_hmm_cl (renamed from hmm_cl), and ZG_hmm_fit
%
% ---INPUTS:
% y, the input time series
% trainp, the proportion of data to train on, 0 < trainp < 1 (default: 0.8)
% numStates, the number of states in the HMM (default: 3)
%
% ---OUTPUTS:
% Mu_1, Mu_2, Mu_3: the state means, sorted from lowest to highest (one per state)
% meanMu, rangeMu, maxMu, minMu: mean, range, maximum, and minimum of the state means
% Cov: the variance of the emissions, shared by all states
% Pmeandiag: mean of the diagonal of the transition matrix (probability of staying
%       in a state)
% stdmeanP: standard deviation across states of the mean probability of moving into
%       each state
% maxP, stdP: maximum and standard deviation of the transition probabilities
% LLtrainpersample: the highest log-likelihood per sample reached on the training part
% LLtestpersample: the log-likelihood per sample of the test part
% LLdifference: LLtestpersample - LLtrainpersample
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

% Check required function files exist:
if ~exist('ZG_hmm', 'file') || ~exist('ZG_hmm_cl', 'file') || ~exist('ZG_hmm_fit', 'file')
	error('Could not find the required HMM fitting functions (Zoubin Gharamani''s code)');
end

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
N = length(y); % number of samples in time series

if nargin < 2 || isempty(trainp)
	fprintf(1, 'Training on 80%% of the data by default\n');
	trainp = 0.8; % train on 80% of the data
end

if nargin < 3 || isempty(numStates)
	fprintf(1, 'Using 3 states by default\n');
	numStates = 3; % use 3 states
end

% ------------------------------------------------------------------------------
%% Train the HMM
% ------------------------------------------------------------------------------
% Divide up dataset into training (ytrain) and test (ytest) portions
Ntrain = floor(trainp * N);
if Ntrain == N
	error('No data for test set for HMM fitting?!')
end

ytrain = y(1:Ntrain);
if Ntrain < N
	ytest = y(Ntrain + 1:end);
	Ntest = length(ytest);
end

% Train HMM with <numStates> states for 30 cycles of EM (or until
% convergence), from fixed starting points; default termination tolerance
[Mu, Cov, P, Pi, LL] = ZG_hmm_fit(ytrain, numStates, 30);

% ------------------------------------------------------------------------------
%% Output statistics on the training
% ------------------------------------------------------------------------------

% Mean vector, Mu
Musort = sort(Mu, 'ascend');
for i = 1:length(Mu)
	out.(sprintf('Mu_%u', i)) = Musort(i); % use dynamic field referencing
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
out.stdP = std(P(:));

% Within-sample log-likelihood
out.LLtrainpersample = max(LL) / Ntrain; % loglikelihood per sample

% ------------------------------------------------------------------------------
%% Calculate log likelihood for the test data
% ------------------------------------------------------------------------------
lik = ZG_hmm_cl(ytest, Ntest, numStates, Mu, Cov, P, Pi);

out.LLtestpersample = lik / Ntest;

out.LLdifference = out.LLtestpersample - out.LLtrainpersample;

end
