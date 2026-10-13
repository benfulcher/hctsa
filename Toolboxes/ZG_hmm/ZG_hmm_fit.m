function [Mu, Cov, P, Pi, LL] = ZG_hmm_fit(X, K, cyc, rhos, floorFrac)
% Deterministic fit of a Gaussian-emission hidden Markov model to a scalar time series
%
% Runs the Baum-Welch fit of ZG_hmm from six fixed starting points and returns the
% fit with the highest final log likelihood on the data X. The starts have a shared
% variance equal to the variance of X, equal initial-state probabilities, and a
% transition matrix with probability rho of staying in a state (and (1-rho)/(K-1) of
% moving to each other one), for each rho in rhos, and two placements of the K state
% means: (i) at the (k-1/2)/K quantiles of X (X sorted, element ceil(N*(k-1/2)/K)),
% and (ii) evenly spaced from mean(X) - std(X) to mean(X) + std(X). The variance
% (shared by all states) is not allowed to fall below floorFrac times the variance of X.
%
% X - N x 1 data
% K - number of states
% cyc - maximum number of cycles of Baum-Welch (default 30)
% rhos - the probabilities of staying in a state (default [0.9, 0.5, 0.99])
% floorFrac - the variance floor, as a proportion of the variance of X (default 0.01)
%
% Mu, Cov, P, Pi, LL - as for ZG_hmm
%
% Reproducible: no random numbers are used.

% ------------------------------------------------------------------------------
if nargin < 3 || isempty(cyc), cyc = 30; end
if nargin < 4 || isempty(rhos), rhos = [0.9, 0.5, 0.99]; end
if nargin < 5 || isempty(floorFrac), floorFrac = 0.01; end

X = X(:);
N = length(X);
v = var(X);
xs = sort(X);
meanSets = [xs(ceil(N * ((1:K)' - 0.5) / K)), mean(X) + std(X) * linspace(-1, 1, K)'];
init.Cov = v;
init.Pi = ones(1, K) / K;

bestLL = -Inf;
Mu = NaN(K, 1); Cov = NaN; P = NaN(K); Pi = NaN(1, K); LL = NaN;
for r = 1:2 * length(rhos)
	init.Mu = meanSets(:, 1 + (r > length(rhos)));
	rho = rhos(1 + mod(r - 1, length(rhos)));
	if K > 1
		init.P = (1 - rho) / (K - 1) * ones(K) + (rho - (1 - rho) / (K - 1)) * eye(K);
	else
		init.P = 1;
	end
	[Mu_r, Cov_r, P_r, Pi_r, LL_r] = ZG_hmm(X, N, K, cyc, [], init, floorFrac * v);
	% (a fit with a non-finite log likelihood or parameters is discarded)
	if all(isfinite([Mu_r(:); Cov_r; P_r(:); Pi_r(:); LL_r(end)])) && LL_r(end) > bestLL
		bestLL = LL_r(end);
		Mu = Mu_r; Cov = Cov_r; P = P_r; Pi = Pi_r; LL = LL_r;
	end
end

end
