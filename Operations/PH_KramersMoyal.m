function out = PH_KramersMoyal(y, tau, numBins)
% PH_KramersMoyal   State-dependent drift and diffusion (Kramers-Moyal coefficients).
%
% Treats the time series as a sampled Langevin process,
%   dx = D1(x) dt + sqrt(2 D2(x)) dW,
% and estimates the conditional-moment (Kramers-Moyal) coefficients as a
% function of the current level x, from increments over tau samples:
%   D1(x) = E[dx | x] / tau             (drift: the deterministic restoring force)
%   D2(x) = E[dx^2 | x] / (2 tau)       (diffusion: the local noise intensity)
%   D4(x) = E[dx^4 | x] / (24 tau)
% with conditioning on x done in equiprobable (quantile) bins.
%
% The features summarize the SHAPE of D1 and D2 across levels:
% - a linear D1 with a constant D2 is an Ornstein-Uhlenbeck (AR(1)-like) process;
% - a cubic D1 (restoring force growing faster than linearly, or several fixed
%   points) indicates a nonlinear potential, e.g., bistability;
% - a D2 that varies with x indicates multiplicative (state-dependent) noise on
%   the scale the data are measured on (note that a monotone transformation of an
%   Ornstein-Uhlenbeck process also has state-dependent D2, so this is partly a
%   property of the marginal distribution);
% - the Pawula ratio D4/D2^2 is small for continuous diffusion and large when
%   increments are dominated by jumps.
% The linear drift slope is closely related to the lag-tau autocorrelation and is
% reported for reference.
%
% ---INPUTS:
% y, the input time series (z-scored in hctsa)
% tau, the increment lag, in samples; or 'ac' for the first zero-crossing of the
%       autocorrelation function
% numBins, the number of equiprobable bins used to condition on x
%
% ---OUTPUTS:
% driftLin, driftQuad, driftCubic: coefficients of a count-weighted cubic fit of
%       D1(x) (driftCubic < 0 for a stiffening restoring force)
% driftNonlinGain: fraction of the residual variance of a linear drift fit that
%       is removed by the cubic fit (0 = linear drift)
% diffSlope, diffCurv: linear and quadratic coefficients of a quadratic fit of
%       D2(x), each relative to the fitted D2 at x = 0 (0 = additive noise)
% diffTailRatio: mean D2 in the outer bins / mean D2 in the central bins
% pawula: mean over bins of D4/D2^2 (0.5 for Gaussian increments; larger for
%       jumps). Unlike the shape of D2, this is not changed by a monotone
%       transformation of the data, which maps a diffusion to another diffusion.

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
% Check inputs and set defaults:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(tau)
    tau = 1;
end
if nargin < 3 || isempty(numBins)
    numBins = 15;
end
y = y(:);

if ischar(tau)
    switch tau
    case 'ac'
        tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
    otherwise
        error('Unknown time delay ''%s''', tau);
    end
    if isnan(tau)
        out = NaN; return % data-dependent: no correlation length could be estimated
    end
end

x = y(1:end-tau);
dx = y(1+tau:end) - x;
if numel(x) < 20*numBins
    out = NaN; return % data-dependent: too few samples per bin
end

% ------------------------------------------------------------------------------
% Conditional moments in equiprobable bins of x:
% ------------------------------------------------------------------------------
edges = quantile(x, linspace(0, 1, numBins + 1));
edges(1) = -Inf; edges(end) = Inf;
bin = discretize(x, edges);
if any(isnan(bin)) || numel(unique(bin)) < numBins
    out = NaN; return % data-dependent: heavily tied values (quantile edges coincide)
end
n = accumarray(bin, 1, [numBins, 1]);
xc = accumarray(bin, x, [numBins, 1], @median);
D1 = accumarray(bin, dx, [numBins, 1], @mean) / tau;
D2 = accumarray(bin, dx.^2, [numBins, 1], @mean) / (2*tau);
D4 = accumarray(bin, dx.^4, [numBins, 1], @mean) / (24*tau);

% ------------------------------------------------------------------------------
% Drift: linear vs cubic fit (weighted by bin counts, which are ~equal)
% ------------------------------------------------------------------------------
w = sqrt(n);
V3 = [ones(numBins,1), xc, xc.^2, xc.^3];
c3 = (V3 .* w) \ (D1 .* w);
V1 = V3(:, 1:2);
c1 = (V1 .* w) \ (D1 .* w);
res1 = sum(n .* (D1 - V1*c1).^2);
res3 = sum(n .* (D1 - V3*c3).^2);
out.driftLin = c3(2);
out.driftQuad = c3(3);
out.driftCubic = c3(4);
out.driftNonlinGain = 1 - res3/res1;

% ------------------------------------------------------------------------------
% Diffusion: quadratic fit relative to its value at x = 0
% ------------------------------------------------------------------------------
V2 = V3(:, 1:3);
b = (V2 .* w) \ (D2 .* w);
out.diffSlope = b(2)/b(1);
out.diffCurv = b(3)/b(1);
numOuter = max(1, floor(numBins/5));
isOuter = [true(numOuter,1); false(numBins - 2*numOuter, 1); true(numOuter,1)];
isCentral = false(numBins,1);
isCentral(ceil(numBins/2) + (-floor(numOuter/2):floor(numOuter/2))) = true;
out.diffTailRatio = mean(D2(isOuter))/mean(D2(isCentral));

% ------------------------------------------------------------------------------
% Pawula ratio (jumps vs continuous diffusion)
% ------------------------------------------------------------------------------
% (the mean over bins is more reliable than the median: test-retest 0.87 vs 0.76)
out.pawula = mean(D4 ./ D2.^2);

end
