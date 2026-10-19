function varargout = BF_GARCHFit(action, varargin)
% BF_GARCHFit   Maximum-likelihood estimation of GARCH, GJR-GARCH and t-GARCH(P,Q) models.
%
% Specified here from the published model definitions (Bollerslev's GARCH, the threshold
% extension of Glosten, Jagannathan and Runkle, and Bollerslev's Student-t innovations),
% not taken from an estimation toolbox. Used by MF_GARCHfit and MF_GARCHcompare, and
% mirrored in pyhctsa.
%
% Model. For a zero-mean series y_1,...,y_N the conditional variance is
%
%   h_t = w + sum_{i=1}^Q (a_i + g_i*[y_{t-i} < 0]) y_{t-i}^2 + sum_{j=1}^P b_j h_{t-j}
%
% (g_i = 0 for a plain GARCH). The innovations y_t/sqrt(h_t) are standard normal or
% Student-t with nu = 1/kappa degrees of freedom, scaled to unit variance.
%
% Presample. y_s^2 = h_s = B for s <= 0, where B is the mean of y^2 (the variance the
% model would have without clustering); the presample leverage indicator is 1/2.
%
% Parameters and bounds. theta = (w, b_1..b_P, a_1..a_Q, g_1..g_Q, kappa), the last two
% groups only for GJR and t. The feasible set is a polytope:
%   1e-6 <= w <= 10,  b_j >= 0,  a_i >= 0,  a_i + g_i >= 0,
%   sum(a) + sum(g)/2 + sum(b) <= 0.9999 (persistence below one),
%   0.01 <= kappa <= 1/2.2 (2.2 <= nu <= 100).
% The leverage g_i may be negative (down to -a_i), so its sign is informative.
%
% Estimation. A deterministic grid of starting points (total persistence, the share of it
% in the ARCH terms and, for t, the tail weight) is scored by the likelihood. The best 12 grid
% points are taken through a coarse sequence of log-barrier Newton problems (barrier weights
% 1e-3 to 1e-5); the two best are then taken to barrier weight 1e-10 (the Newton polish),
% and the better is the fit. Newton steps use the analytic gradient and a central-difference
% Hessian of it, and flip the sign of negative curvature. Nothing is random.
%
% ---USAGE (the first input selects the action):
% D = BF_GARCHFit('data', y)
%       prepares the (whitened, z-scored) series y
% S = BF_GARCHFit('spec', D, P, Q, isGJR, isT)
%       the model structure (parameter positions, constraints and lagged squares)
% fit = BF_GARCHFit('mle', D, P, Q, isGJR, isT, extraStarts, useGrid)
%       the maximum-likelihood fit; extraStarts is a matrix of starting parameter vectors
%       (columns), useGrid is whether the grid of starts is also used (default: true; the
%       grid is used anyway if there are no extra starts). The fit has fields theta, llf
%       (total log-likelihood), h (conditional variances), spec, dec (the last Newton
%       decrement, small when converged) and k (the number of parameters)
% fit = BF_GARCHFit('null', D, isT)
%       the constant-variance special case (Gaussian in closed form; t by maximum likelihood)
% theta = BF_GARCHFit('embed', theta, specFrom, specTo)
%       a smaller (P,Q) model's parameters in a larger model (new coefficients at zero)
% starts = BF_GARCHFit('tailstarts', gaussFit, specT)
%       starting points of a t model from a Gaussian fit of the same orders
% [lr, hasStructure] = BF_GARCHFit('structure', alt, null, D, isGJR)
%       the likelihood-ratio test of a variance model against constant variance. For Gaussian
%       innovations the statistic is rescaled by 2/(kurtosis - 1) of the standardized
%       residuals (heavy-tailed innovations inflate the Gaussian statistic); for t it is
%       used as it is. Its null distribution is not chi-squared, because a parameter lies on
%       the boundary and another is not identified under the null (D. W. K. Andrews,
%       Econometrica 69, 683-734 (2001), DOI: 10.1111/1468-0262.00210), so the critical
%       values, 4 (6 for GJR), were calibrated by simulation of white noise: 3% to 6% of white-noise
%       series of 100 to 10000 samples exceed 4 (about 5% exceed 6 for GJR).
% errors = BF_GARCHFit('errors', theta, D, spec)
%       standard errors from the outer product of the per-observation scores
% atFloor = BF_GARCHFit('atfloor', theta)
%       whether the variance constant of a fit sits at its lower bound, 1e-6, within 1%. The
%       Gaussian likelihood of a series with long stretches of (near) zero values has no
%       maximum (the variance can shrink without limit where the series is zero), so a fit
%       at the floor has a likelihood and a variance path that depend on the bound rather
%       than on the data
% tailLR = BF_GARCHFit('taillr')
%       the critical value of the likelihood-ratio test of Student-t against Gaussian
%       innovations: 2.7055, the 5% point of the 50:50 mixture of chi-squared(0) and
%       chi-squared(1) (a Gaussian is the limit of infinite degrees of freedom, on the
%       boundary of the parameter space; Self and Liang, J. Am. Stat. Assoc. 82, 605 (1987),
%       DOI: 10.1080/01621459.1987.10478472)
%
% References:
% T. Bollerslev, "Generalized autoregressive conditional heteroskedasticity", J.
% Econometrics 31, 307 (1986). DOI: 10.1016/0304-4076(86)90063-1
% L. R. Glosten, R. Jagannathan and D. E. Runkle, "On the relation between the expected
% value and the volatility of the nominal excess return on stocks", J. Finance 48, 1779
% (1993). DOI: 10.1111/j.1540-6261.1993.tb05128.x
% T. Bollerslev, "A conditionally heteroskedastic time series model for speculative prices
% and rates of return", Rev. Econ. Stat. 69, 542 (1987). DOI: 10.2307/1925546

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

switch action
case 'data'
    varargout{1} = makeData(varargin{1});
case 'mle'
    varargout{1} = mleFit(varargin{:});
case 'null'
    varargout{1} = nullFit(varargin{:});
case 'embed'
    varargout{1} = embedParams(varargin{:});
case 'tailstarts'
    varargout{1} = tailStarts(varargin{:});
case 'structure'
    [varargout{1}, varargout{2}] = structureTest(varargin{:});
case 'errors'
    varargout{1} = opgErrors(varargin{:});
case 'taillr'
    varargout{1} = 2.7055;
case 'atfloor'
    varargout{1} = varargin{1}(1) <= 1e-6 * (1 + 0.01);
case 'spec'
    varargout{1} = makeSpec(varargin{:});
otherwise
    error('Unknown action ''%s''', action);
end

end

% ------------------------------------------------------------------------------
function D = makeData(y)
y = y(:);
D.y = y;
D.N = length(y);
D.y2 = y.^2;
D.B = mean(D.y2);
D.y2g = D.y2 .* (y < 0);
end

% ------------------------------------------------------------------------------
function S = makeSpec(D, P, Q, gjr, tdist)
% Model structure, its linear constraints A*theta <= c, and the lagged squares
W_LO = 1e-6; W_HI = 10; RHO_MAX = 0.9999; KAPPA_LO = 0.01; KAPPA_HI = 1/2.2;
S.P = P; S.Q = Q; S.gjr = logical(gjr); S.t = logical(tdist);
S.ib = 1 + (1:P);
S.ia = 1 + P + (1:Q);
if gjr
    S.ig = 1 + P + Q + (1:Q);
else
    S.ig = [];
end
if tdist
    S.ik = 1 + P + Q + gjr*Q + 1;
else
    S.ik = [];
end
S.k = 1 + P + Q + gjr*Q + tdist;
rows = zeros(0, S.k); rhs = zeros(0, 1);
    function addRow(idx, vals, b)
        r = zeros(1, S.k); r(idx) = vals;
        rows(end+1, :) = r; %#ok<AGROW>
        rhs(end+1, 1) = b; %#ok<AGROW>
    end
addRow(1, -1, -W_LO);
addRow(1, 1, W_HI);
for j = S.ib, addRow(j, -1, 0); end
for i = S.ia, addRow(i, -1, 0); end
for n = 1:length(S.ig), addRow([S.ia(n), S.ig(n)], [-1, -1], 0); end
if P + Q > 0
    addRow([S.ib, S.ia, S.ig], [ones(1, P + Q), 0.5*ones(1, length(S.ig))], RHO_MAX);
end
if tdist
    addRow(S.ik, -1, -KAPPA_LO);
    addRow(S.ik, 1, KAPPA_HI);
end
S.A = rows; S.c = rhs;
% Lagged squares (presample B) and lagged squares of negative values (presample B/2)
N = D.N; B = D.B;
S.L = zeros(N, Q); S.Lg = zeros(N, Q);
for i = 1:Q
    S.L(1:i, i) = B;
    S.L(i+1:N, i) = D.y2(1:N-i);
    S.Lg(1:i, i) = B/2;
    S.Lg(i+1:N, i) = D.y2g(1:N-i);
end
end

% ------------------------------------------------------------------------------
function [h, dh] = condVariance(theta, D, S, wantGrad)
N = D.N; B = D.B; P = S.P; Q = S.Q;
d = theta(1) * ones(N, 1);
if Q > 0
    d = d + S.L * theta(S.ia);
end
if S.gjr
    d = d + S.Lg * theta(S.ig);
end
if P > 0
    b = theta(S.ib);
    A = [1; -b];
    zi = B * flipud(cumsum(flipud(b)));   % filter state for a presample of B
    h = filter(1, A, d, zi);
else
    A = 1;
    h = d;
end
dh = [];
if ~wantGrad
    return
end
M = zeros(N, S.k);
M(:, 1) = 1;
for j = 1:P
    hl = zeros(N, 1);
    hl(1:j) = B;
    hl(j+1:N) = h(1:N-j);
    M(:, S.ib(j)) = hl;
end
for i = 1:Q
    M(:, S.ia(i)) = S.L(:, i);
end
for i = 1:length(S.ig)
    M(:, S.ig(i)) = S.Lg(:, i);
end
if P > 0
    dh = filter(1, A, M);
else
    dh = M;
end
dh(:, S.k - S.t + 1:S.k) = 0;   % (the kappa column carries no variance derivative)
end

% ------------------------------------------------------------------------------
function [f, g, sc] = nll(theta, D, S, wantGrad, wantScores)
% Negative log-likelihood per observation, its gradient and per-observation scores
y2 = D.y2; N = D.N;
g = []; sc = [];
[h, dh] = condVariance(theta, D, S, wantGrad);
if any(h <= 0) || any(~isfinite(h))
    f = Inf;
    return
end
if ~S.t
    f = 0.5 * (log(2*pi) + mean(log(h) + y2 ./ h));
    if wantGrad
        w1 = 0.5 * (1 ./ h - y2 ./ h.^2);
        g = (w1' * dh)' / N;
        if wantScores
            sc = w1 .* dh;
        end
    end
    return
end
kappa = theta(S.ik); nu = 1 / kappa;
w = y2 ./ (h * (nu - 2));
lw = log1p(w);
c = gammaln((nu + 1)/2) - gammaln(nu/2) - 0.5 * log(pi * (nu - 2));
f = -c + 0.5 * mean(log(h)) + 0.5 * (nu + 1) * mean(lw);
if wantGrad
    w1 = 0.5 ./ h .* (1 - (nu + 1) * w ./ (1 + w));
    g = (w1' * dh)' / N;
    cnu = 0.5 * psi((nu + 1)/2) - 0.5 * psi(nu/2) - 0.5 / (nu - 2);
    dnuObs = -cnu + 0.5 * lw - 0.5 * (nu + 1) * w ./ ((nu - 2) * (1 + w));
    g(S.ik) = mean(dnuObs) * (-nu^2);
    if wantScores
        sc = w1 .* dh;
        sc(:, S.ik) = dnuObs * (-nu^2);
    end
end
end

% ------------------------------------------------------------------------------
function room = slackSteps(theta, S)
% The largest step in each coordinate that keeps every constraint satisfied
s = S.c - S.A * theta;
room = Inf(length(theta), 1);
for r = 1:size(S.A, 1)
    idx = find(S.A(r, :) ~= 0);
    for i = idx
        room(i) = min(room(i), s(r) / abs(S.A(r, i)));
    end
end
end

% ------------------------------------------------------------------------------
function H = numHessian(theta, D, S)
k = length(theta);
room = slackSteps(theta, S);
H = zeros(k, k);
g0 = [];
for i = 1:k
    step = max(min(1e-5 * max(1, abs(theta(i))), 0.4 * room(i)), 1e-13);
    e = zeros(k, 1); e(i) = step;
    [~, gp] = nll(theta + e, D, S, true, false);
    [~, gm] = nll(theta - e, D, S, true, false);
    if isempty(gp) || isempty(gm)
        if isempty(g0)
            [~, g0] = nll(theta, D, S, true, false);
        end
        if isempty(gp), gp = g0; end
        if isempty(gm), gm = g0; end
    end
    H(:, i) = (gp - gm) / (2 * step);
end
H = 0.5 * (H + H');
end

% ------------------------------------------------------------------------------
function [phi, g, s] = barrier(theta, D, S, mu, wantGrad)
% (wantGrad false: the value only, as needed in the line search)
if nargin < 5, wantGrad = true; end
s = S.c - S.A * theta;
g = [];
if any(s <= 0)
    phi = Inf;
    return
end
if ~wantGrad
    f = nll(theta, D, S, false, false);
    phi = f - mu * sum(log(s));
    return
end
[f, gf] = nll(theta, D, S, true, false);
if isempty(gf)
    phi = Inf;
    return
end
phi = f - mu * sum(log(s));
g = gf + mu * (S.A' * (1 ./ s));
end

% ------------------------------------------------------------------------------
function [theta, dec] = newtonBarrier(theta, D, S, mus)
maxIt = 60; tol = 1e-13;
dec = Inf;
for mu = mus
    for it = 1:maxIt
        [phi, g, s] = barrier(theta, D, S, mu);
        H = numHessian(theta, D, S) + mu * (S.A' * diag(1 ./ s.^2) * S.A);
        H = 0.5 * (H + H');
        [V, E] = eig(H);
        ev = diag(E);
        ev = max(abs(ev), 1e-9 * max(1, max(abs(ev))));
        step = -V * ((V' * g) ./ ev);
        dec = -g' * step;
        if dec < tol
            break
        end
        Ad = S.A * step;
        pos = Ad > 0;
        if ~any(pos)
            a = 1;
        else
            a = min(1, 0.995 * min(s(pos) ./ Ad(pos)));
        end
        ok = false;
        for ls = 1:60
            if barrier(theta + a * step, D, S, mu, false) <= phi - 1e-4 * a * dec
                ok = true;
                break
            end
            a = a * 0.5;
        end
        if ~ok
            break
        end
        theta = theta + a * step;
    end
end
end

% ------------------------------------------------------------------------------
function pts = startGrid(D, S)
% Grid of starting points (columns), sorted by negative log-likelihood (ties in grid order)
P = S.P; Q = S.Q; B = D.B;
if P + Q == 0
    rhos = 0;
else
    rhos = [0.3, 0.6, 0.9, 0.98];
end
if P > 0 && Q > 0
    shares = [0.1, 0.25, 0.5];
else
    shares = 1;
end
if S.t
    kappas = [0.02, 0.1, 0.2];
else
    kappas = NaN;
end
pts = zeros(S.k, 0);
for rho = rhos
    for sh = shares
        th = zeros(S.k, 1);
        if P > 0 && Q > 0
            arch = rho * sh; gar = rho * (1 - sh);
        elseif Q > 0
            arch = rho; gar = 0;
        else
            arch = 0; gar = rho;
        end
        if P > 0, th(S.ib) = gar / P; end
        if Q > 0, th(S.ia) = arch / Q; end
        th(1) = max(B * (1 - rho), 2e-6);
        for kp = kappas
            t2 = th;
            if ~isnan(kp), t2(S.ik) = kp; end
            pts(:, end+1) = t2; %#ok<AGROW>
        end
    end
end
sc = zeros(1, size(pts, 2));
for n = 1:size(pts, 2)
    sc(n) = nll(pts(:, n), D, S, false, false);
end
[~, order] = sort(sc);
pts = pts(:, order);
end

% ------------------------------------------------------------------------------
function th = makeInterior(theta, S)
% Move a start strictly inside the constraints (towards a fixed interior reference point)
ref = zeros(S.k, 1);
ref(1) = 0.5;
ref(S.ib) = 0.05 / max(S.P, 1);
ref(S.ia) = 0.05 / max(S.Q, 1);
if S.t, ref(S.ik) = 0.1; end
th = theta;
for n = 1:50
    if all(S.c - S.A * th > 5e-7)
        return
    end
    th = 0.9 * th + 0.1 * ref;
end
end

% ------------------------------------------------------------------------------
function fit = mleFit(D, P, Q, gjr, tdist, extraStarts, useGrid, npolish)
if nargin < 6, extraStarts = zeros(0, 0); end
if nargin < 7 || isempty(useGrid), useGrid = true; end
if nargin < 8, npolish = 2; end
muCoarse = [1e-3, 1e-4, 1e-5];
muFine = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10];
S = makeSpec(D, P, Q, gjr, tdist);
starts = extraStarts;
if useGrid || isempty(extraStarts)
    grid = startGrid(D, S);
    starts = [starts, grid(:, 1:min(12, size(grid, 2)))]; % (the best 12 by likelihood)
end
nS = size(starts, 2);
phis = zeros(1, nS); ths = zeros(S.k, nS);
for n = 1:nS
    th = newtonBarrier(makeInterior(starts(:, n), S), D, S, muCoarse);
    ths(:, n) = th;
    phis(n) = barrier(th, D, S, muCoarse(end));
end
[~, order] = sort(phis);
bestF = Inf; bestTh = []; bestDec = NaN;
for n = order(1:min(npolish, nS))
    [th, dec] = newtonBarrier(ths(:, n), D, S, muFine);
    f = nll(th, D, S, false, false);
    if f < bestF
        bestF = f; bestTh = th; bestDec = dec;
    end
end
fit.theta = bestTh;
fit.llf = -bestF * D.N;
fit.h = condVariance(bestTh, D, S, false);
fit.spec = S;
fit.dec = bestDec;
fit.k = S.k;
end

% ------------------------------------------------------------------------------
function fit = nullFit(D, tdist)
if ~tdist
    fit.theta = D.B;
    fit.llf = -0.5 * D.N * (log(2*pi) + log(D.B) + 1);
    fit.h = D.B * ones(D.N, 1);
    fit.spec = makeSpec(D, 0, 0, false, false);
    fit.dec = 0;
    fit.k = 1;
else
    fit = mleFit(D, 0, 0, false, true);
end
end

% ------------------------------------------------------------------------------
function th = embedParams(theta, sFrom, sTo)
th = zeros(sTo.k, 1);
th(1) = theta(1);
th(sTo.ib(1:sFrom.P)) = theta(sFrom.ib);
th(sTo.ia(1:sFrom.Q)) = theta(sFrom.ia);
if ~isempty(sFrom.ig)
    th(sTo.ig(1:length(sFrom.ig))) = theta(sFrom.ig);
end
if sFrom.t && sTo.t
    th(sTo.ik) = theta(sFrom.ik);
end
end

% ------------------------------------------------------------------------------
function starts = tailStarts(gaussFit, sT)
starts = zeros(sT.k, 3);
kps = [0.02, 0.1, 0.2];
for n = 1:3
    th = embedParams(gaussFit.theta, gaussFit.spec, sT);
    th(sT.ik) = kps(n);
    starts(:, n) = th;
end
end

% ------------------------------------------------------------------------------
function [lr, has] = structureTest(alt, null, D, gjr)
if nargin < 4, gjr = false; end
lr = max(2 * (alt.llf - null.llf), 0);
if ~alt.spec.t
    z2 = D.y2 ./ alt.h;
    lr = lr * 2 / max(mean(z2.^2) / mean(z2)^2 - 1, 1e-3);
end
if gjr
    crit = 6;
else
    crit = 4;
end
has = lr >= crit;
end

% ------------------------------------------------------------------------------
function errors = opgErrors(theta, D, S)
[~, ~, sc] = nll(theta, D, S, true, true);
if isempty(sc)
    errors = NaN(length(theta), 1);
    return
end
errors = sqrt(abs(diag(BF_Pinv(sc' * sc, 1e-10))));
end
