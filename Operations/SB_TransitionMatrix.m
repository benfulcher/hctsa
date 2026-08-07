function out = SB_TransitionMatrix(y, howtocg, numGroups, tau)
% SB_TransitionMatrix  Transition probabilities between time-series states.
%
% The time series is coarse-grained according to a given method.
%
% The input time series is transformed into a symbolic string using an
% equiprobable alphabet of numGroups letters. The transition probabilities are
% calculated at a lag tau.
%
% Related to the idea of quantile graphs from time series.
% cf. Andriana et al. (2011). Duality between Time Series and Networks. PLoS ONE.
% https://doi.org/10.1371/journal.pone.0023378
%
% ---INPUTS:
% y, the input time series
%
% howtocg, the method of discretization: 'quantile' (equiprobable, the
%           default) or 'updown' (a true binary up/down split by the sign of
%           each increment -- NOT equiprobable, and requires numGroups=2; see
%           SB_CoarseGrain.m). Other SB_CoarseGrain methods could be
%           incorporated in future.
%
% numGroups: number of groups in the course-graining
%
% tau: analyze transition matricies corresponding to this lag. We
%      could either downsample the time series at this lag and then do the
%      discretization as normal, or do the discretization and then just
%      look at this dicrete lag. Here we do the former. Can also set tau to 'ac'
%      to set tau to the first zero-crossing of the autocorrelation function.
%
% ---OUTPUTS: include the transition probabilities themselves, as well as the trace
% of the transition matrix, measures of asymmetry, and eigenvalues of the
% transition matrix.

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
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(howtocg)
	howtocg = 'quantile';
end
if nargin < 3 || isempty(numGroups)
	numGroups = 2;
end
if numGroups < 2
	error('Too few groups for coarse-graining')
end
if nargin < 4 || isempty(tau)
	tau = 1;
end
if strcmp(tau, 'ac') % determine tau from first zero of autocorrelation
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
end
if isnan(tau)
	error('Time series too short to estimate tau');
end

if tau > 1 % calculate transition matrix at a non-unit lag
	% downsample at rate 1:tau
	y = resample(y, 1, tau);
end

% ------------------------------------------------------------------------------
%% (((1))) Discretize the time series to a symbolic string
% ------------------------------------------------------------------------------
yth = SB_CoarseGrain(y, howtocg, numGroups);

% At this point we should have:
% (*) yth: a thresholded y containing integers from 1 to numGroups

if size(yth, 2) > size(yth, 1)
	yth = yth';
end

% N is taken AFTER coarse-graining, not length(y): some coarse-graining
% methods (e.g. 'updown', which internally differences the series) return a
% yth shorter than y, and using length(y) here would silently mis-normalize
% T below.
N = length(yth);

% ------------------------------------------------------------------------------
%% (((2))) Compute the tau-step transition matrix
%               (Markov for tau = 1)
% ------------------------------------------------------------------------------
% Probably implemented already, but I'll do it myself
T = zeros(numGroups); % probability of transition from state i -> state j
for i = 1:numGroups
	ri = (yth == i); % indices where the time series is in state i
	if sum(ri) == 0 % is never in state i
		T(i, :) = 0; % all transition probabilities are zero (could be NaN)
	else
		% Indices of states immediately following a state i:
		ri_next = [false; ri(1:end - 1)];
		% Compute transitions from state i to each of the states j:
		for j = 1:numGroups
			T(i, j) = sum(yth(ri_next) == j); % the next element is of this class
		end
	end
end

% Normalize from counts to probabilities:
T = T / (N - 1); % N-1 is appropriate because it's a 1-time transition matrix

% ------------------------------------------------------------------------------
%% (((3))) Output measures from the transition matrix
% ------------------------------------------------------------------------------
% (i) Raw values of the transition matrix
% [this has to be done bulkily (only for numGroups = 2,3)]:
if numGroups == 2 % return all elements of T
	for i = 1:4
		out.(sprintf('T%u', i)) = T(i);
	end
elseif numGroups == 3 % return all elements of T
	for i = 1:9
		out.(sprintf('T%u', i)) = T(i);
	end
elseif numGroups > 3 % return just diagonal elements of T
	for i = 1:numGroups
		out.(sprintf('TD%u', i)) = T(i, i);
	end
end

% (ii) Measures on the diagonal
out.ondiag = sum(diag(T)); % trace
out.stddiag = std(diag(T)); % std of diagonal elements

% (iii) Measures of symmetry:
out.symdiff = sum(sum(abs((T - T')))); % sum of differences of individual elements
out.symsumdiff = sum(sum(tril(T, -1))) - sum(sum(triu(T, +1))); % difference in sums of upper and lower
% triangular parts of T

% Kullback-Leibler divergence between T and its transpose T': a
% reversal-asymmetry ("irreversibility") measure with a direct
% information-theoretic interpretation (related to entropy production rate
% for a Markov chain), zero iff T satisfies detailed balance (T(i,j)=T(j,i)
% for all i,j, i.e. forward and backward transition rates match exactly).
% Restricted to pairs where both T(i,j) and T(j,i) are nonzero (a one-sided
% pair, observed one direction but never the other, would otherwise give an
% infinite contribution) -- this keeps the measure finite and well-behaved
% for the small numGroups used here, at the cost of (deliberately)
% understating asymmetry from never-observed reverse transitions.
Tt = T'; % Tt(i,j) = T(j,i), the reverse-direction transition probability
klMask = (T > 0) & (Tt > 0);
out.transKLdiv = sum(T(klMask) .* log(T(klMask) ./ Tt(klMask)));

% (iv) Measures from eigenvalues of T
eigT = eig(T);
out.stdeig = std(eigT); % std of eigenvalues
out.maxeig = max(real(eigT)); % maximum eigenvalue
out.mineig = min(real(eigT)); % minimum eigenvalue
% mean eigenvalue is equivalent to trace
% (ought to be always zero? Not necessary to measure:)
out.maximeig = max(imag(eigT)); % maximum imaginary part of eigenvalues

% Second-largest eigenvalue and spectral gap: for the 'quantile' coarse-
% graining (equiprobable by construction), the marginal (row/column sums) is
% close to uniform, which pins the leading eigenvalue near 1/numGroups
% regardless of temporal structure -- it's the second-largest eigenvalue that
% reflects the Markov chain's mixing rate (a smaller gap indicates slower
% relaxation, i.e. longer memory). This doesn't hold for SB_CoarseGrain's
% non-equiprobable methods ('updown', 'embed2quadrants'/'embed2octants'),
% where the marginal (and hence the leading eigenvalue) can vary meaningfully
% with the data -- these fields aren't registered for the numGroups=2 case
% regardless of coarse-graining method (see the mops file: for 'quantile' it
% would duplicate SB_MotifTwo_median, for 'updown' SB_MotifTwo_diff).
realEig = sort(real(eigT), 'descend');
if numGroups >= 2
	out.secondeig = realEig(2);
	out.specgap = out.maxeig - out.secondeig;
else
	out.secondeig = NaN;
	out.specgap = NaN;
end

% (vii) Transition (conditional) entropy rate, H(X_{t+1}|X_t), treating T as a
% first-order Markov approximation: H(joint) - H(marginal)
rowSums = sum(T, 2); % marginal distribution over states at time t
Hjoint = -sum(T(T > 0) .* log(T(T > 0)));
Hmarginal = -sum(rowSums(rowSums > 0) .* log(rowSums(rowSums > 0)));
% Miller-Madow correction. Each plug-in entropy is biased low by (M-1)/(2n) for
% M occupied bins, so their difference carries a residual bias of
% -(Mjoint - Mmarginal)/(2n) -- a pure function of the number of transitions,
% and hence of time-series length. Adding it back removes the leading 1/n term.
numTransitions = N - 1; % T was normalized by N-1 above
Mjoint = sum(T(:) > 0);
Mmarginal = sum(rowSums > 0);
out.transEntropy = Hjoint - Hmarginal + (Mjoint - Mmarginal) / (2 * numTransitions);

% -------------------------------------------------------------------------------
% (v) Measures from covariance matrix:
covT = cov(T);
out.sumdiagcov = trace(covT); % trace of covariance matrix
% This is equivalent to the sum of column variances: sum([var(T(:,1)),var(T(:,2)),var(T(:,3))])
% (or, similarly, to the sum of row variances): sum([var(T(1,:)),var(T(2,:)),var(T(3,:))])

% (vi) Eigenvalues of covariance matrix
% (mean eigenvalue of covariance matrix equivalent to trace of covariance matrix).
% (these measures don't make much sense in the case of 2 groups):
eigcovT = eig(covT);
out.stdeigcov = std(eigcovT); % std of eigenvalues of covariance matrix
out.maxeigcov = max(eigcovT); % max eigenvalue of covariance matrix
out.mineigcov = min(eigcovT); % min eigenvalue of covariance matrix

end
