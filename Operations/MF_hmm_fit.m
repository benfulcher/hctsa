% MF_hmm_fit
% 
% Uses Zoubin Gharamani's implementation of HMMs for real-valued Gaussian
% observations:
% http://www.gatsby.ucl.ac.uk/~zoubin/software.html
% or, specifically:
% http://www.gatsby.ucl.ac.uk/~zoubin/software/hmm.tar.gz
% 
% Uses ZG_hmm (renamed from hmm) and ZG_hmm_cl (renamed from hmm_cl)
% 
% INPUTS:
% y, the input time series
% trainp, the proportion of data to train on, 0 < trainp < 1
% nstates, the number of states in the HMM
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = MF_hmm_fit(y,trainp,nstates)
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

% Mean vector, Mu
Musort = sort(Mu,'ascend');
for i = 1:length(Mu)
    eval(sprintf('out.Mu_%u = Musort(%u);',i,i));
end
out.meanMu = mean(Mu);
out.rangeMu = max(Mu) - min(Mu);
out.maxMu = max(Mu);
out.minMu = min(Mu);

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