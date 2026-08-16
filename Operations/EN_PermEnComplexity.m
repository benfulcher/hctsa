function out = EN_PermEnComplexity(y, m, tau)
% EN_PermEnComplexity   Jensen-Shannon statistical complexity of ordinal patterns.
%
% Computes the Bandt-Pompe ordinal-pattern distribution (as in EN_PermEn) and
% pairs its normalized Shannon entropy with the Jensen-Shannon statistical
% complexity of Rosso, Larrondo, Martin, Plastino & Fuentes, "Distinguishing
% Noise from Chaos", Phys. Rev. Lett. 99, 154102 (2007) -- the
% entropy-complexity plane used to separate chaotic, stochastic and periodic
% dynamics that can look alike under entropy alone.
%
% Entropy is near its extremes (0 or log(m!)) for both fully ordered *and*
% fully random sequences. The statistical complexity
% C = Q_J[P,P_uniform] . H[P] is instead close to zero at both those extremes
% and peaks for structured-but-disordered ('chaotic') ordinal-pattern
% distributions -- a distinct axis of information from entropy alone, not
% captured elsewhere in hctsa. cf. EN_PermEn for the plain Bandt-Pompe
% permutation entropy this complements (and whose normalized entropy,
% normPermEn, this function's hNorm reproduces).
%
% cf. Martin, Plastino & Rosso, Physica A 369(2) 439 (2006) for the Q_0
% normalization; Lamberti, Martin, Plastino & Rosso, Physica A 334(1-2) 119
% (2004) for the Jensen-Shannon statistical complexity construction.
%
% Note: at m = 2 there are only two ordinal states, so H and C are both
% unimodal, symmetric functions of a single probability -- they are then
% forced to be near-perfect reparameterizations of one another (|r| ~ 1)
% regardless of the input data, making jsComplexity redundant with plain
% permutation entropy at that order. m = 3 was also found redundant with
% existing normPermEn/motif fields on real-world data (Empirical1000,
% r up to 0.97); only m = 4 and m = 5 are registered in the default feature
% set.
%
% ---INPUTS:
% y, the input time series
% m, the embedding dimension (order of the ordinal patterns)
% tau, the time-delay for the embedding
%
% ---OUTPUTS:
% hNorm, normalized Shannon entropy of the ordinal-pattern distribution,
%        H[P] = S[P]/log2(m!), in [0,1]
% jsComplexity, the Jensen-Shannon statistical complexity,
%               C[P] = Q_J[P,P_uniform] . H[P], in [0,1]

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
%% Check inputs & set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(m)
	m = 2; % order 2
end
if nargin < 3 || isempty(tau)
	tau = 1;
end

% ------------------------------------------------------------------------------
%% Embed the signal and count ordinal patterns (as in EN_PermEn):
% ------------------------------------------------------------------------------
x = BF_Embed(y, tau, m, false);
Nx = size(x, 1); % number of embedding vectors produced
if Nx < 5 % need at least 5 embedding vectors to actually do a computation
	% Data-dependent (series too short for the requested tau/m), not a code error.
	warning('Time series too short to embed');
	out.hNorm = NaN;
	out.jsComplexity = NaN;
	return
end

numPerms = factorial(m); % = m!
permIdx = BF_OrdinalPatternRank(x); % index in 1:m! for each embedding vector; cf. EN_PermEn.m
countPerms = accumarray(permIdx, 1, [numPerms, 1]);

p = countPerms / Nx; % ordinal-pattern probability distribution, P

% ------------------------------------------------------------------------------
%% Normalized Shannon entropy of P:
% ------------------------------------------------------------------------------
p_0 = p(p > 0); % makes log(0) = 0
S_P = -sum(p_0 .* log2(p_0));
S_max = log2(numPerms); % entropy of the uniform distribution, Pe
out.hNorm = S_P / S_max;

% ------------------------------------------------------------------------------
%% Jensen-Shannon statistical complexity of P relative to Pe (Rosso et al., 2007):
% ------------------------------------------------------------------------------
pe = 1 / numPerms; % uniform reference distribution, Pe (every entry equal)
pMix = (p + pe) / 2; % (P + Pe)/2 -- every entry > 0 since pe > 0
S_mix = -sum(pMix .* log2(pMix));
S_Pe = S_max; % entropy of Pe is, by construction, log2(numPerms)

% Jensen-Shannon divergence between P and Pe:
JSdiv = S_mix - S_P / 2 - S_Pe / 2;

% Normalizing constant so that Q_J in [0,1] (Martin, Plastino & Rosso, 2006),
% attained in the limit of P a point mass (single ordinal pattern):
N = numPerms;
Q0 = -2 / (((N + 1) / N) * log2(N + 1) - 2 * log2(2 * N) + log2(N));

QJ = Q0 * JSdiv;
out.jsComplexity = QJ * out.hNorm;

end
