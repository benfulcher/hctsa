function m = BF_SubspaceID(y, order)
% BF_SubspaceID   Subspace (N4SID) identification of a state-space model of a time series.
%
% Identifies the discrete-time innovations-form model
%   x(t+1) = A x(t) + K e(t),   y(t) = C x(t) + e(t)
% of a scalar series by the N4SID subspace algorithm with (regularized) canonical-variate
% weighting: orthogonal projections of the future outputs onto the past; a singular value
% decomposition weighted by (Yf*Yf' + rho*I)^(-1/2), with Yf the block of future outputs
% and rho = 1e-3 times the mean eigenvalue of Yf*Yf', for the extended observability
% matrix; A and C from its shift structure; the noise covariances from the residuals of
% the state regression; and the Kalman gain K from the filter Riccati equation. A is
% stabilized first (eigenvalues outside the unit circle are reflected to 1/conj(lambda)).
% The sign of each state is fixed so that C has no negative entries. The loss is the mean
% squared one-step prediction error over the samples predicted from a diffuse initial
% state (BF_SSPredict; the first n are not).
%
% The past and future horizons are equal, s = r = max(n + 1, p), with p the AIC order of
% an AR model of the series (BF_ARorderAIC), reduced (not below n + 1) until the data
% matrices have at least as many columns as rows, N - 2s + 1 >= 2s + 1.
%
%---INPUTS:
% y, the time series
% order, the model order n, or 'best': the order from 1 to 10 minimizing
%       Akaike's information criterion log(max(loss, 1e-10*mean(y.^2))) + 2*(2n)/N,
%       all fitted with the horizons for n = 10 (orders whose fit fails are skipped)
%
%---OUTPUT:
% m, a structure with fields A (n x n), K (n x 1), C (1 x n), loss, and order.
%
% Errors if the series is too short, not finite, or constant, or if the Riccati
% equation has no stabilizing solution.
%
%---NOTES:
% P. Van Overschee and B. De Moor, "N4SID: Subspace algorithms for the identification of
% combined deterministic-stochastic systems", Automatica 30, 75 (1994), and Subspace
% Identification for Linear Systems, Kluwer (1996); W. E. Larimore, "Canonical variate
% analysis in identification, filtering and adaptive control", Proc. 29th IEEE CDC, 596
% (1990); L. Ljung, System Identification: Theory for the User, 2nd ed., Prentice Hall
% (1999), sec. 10.6. Riccati equation: structured doubling algorithm, E. K.-W. Chu,
% H.-Y. Fan and W.-W. Lin, "A structure-preserving doubling algorithm for continuous-time
% algebraic Riccati equations", Linear Algebra Appl. 396, 55 (2005), in its discrete-time
% form (W.-W. Lin and S.-F. Xu, SIAM J. Matrix Anal. Appl. 28, 26 (2006)). Pseudoinverses
% use BF_Pinv (relative singular-value cutoff 1e-10).

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

y = y(:);
N = length(y);
if any(~isfinite(y))
    error('The series must be finite');
end
if range(y) == 0
    error('The series is constant');
end

if ischar(order)
    if ~strcmp(order, 'best')
        error('Unknown order ''%s''', order);
    end
    s = SUB_horizon(y, 10);
    best = []; bestCrit = Inf; lastErr = [];
    floorL = 1e-10 * mean(y.^2);
    for n = 1:10
        try
            fit = SUB_core(y, n, s, s);
        catch lastErr
            continue
        end
        crit = log(max(fit.loss, floorL)) + 2 * 2 * n / N;
        if crit < bestCrit
            best = fit; bestCrit = crit;
        end
    end
    if isempty(best)
        if ~isempty(lastErr)
            rethrow(lastErr);
        end
        error('The series is too short to fit a state-space model');
    end
    m = best;
else
    s = SUB_horizon(y, order);
    m = SUB_core(y, order, s, s);
end

end
% ------------------------------------------------------------------------------
function s = SUB_horizon(y, order)
    % Past (= future) horizon from the AIC order of an AR model
    N = length(y);
    s = max(order + 1, BF_ARorderAIC(y));
    while s > order + 1 && N - 2*s + 1 < 2*s + 1
        s = s - 1;
    end
end
% ------------------------------------------------------------------------------
function fit = SUB_core(y, n, r, s)
    N = length(y);
    j = N - r - s + 1;
    if n < 1 || s < n || r < n + 1 || j < r + s + 1
        error('The series is too short for the requested order and horizons');
    end
    tt = s + (1:j); % indices of the first future sample of each column
    Ypast = y(tt - 1 - (0:s-1)'); % s x j: y(t-1), ..., y(t-s)
    Yfut = y(tt + (0:r-1)'); % r x j: y(t), ..., y(t+r-1)
    YpastExt = [Yfut(1, :); Ypast];
    YfutMinus = Yfut(2:end, :);

    O = (Yfut * BF_Pinv(Ypast)) * Ypast;
    Ominus = (YfutMinus * BF_Pinv(YpastExt)) * YpastExt;

    % Canonical-variate weighting by the inverse square root of the future covariance,
    % regularized by a ridge of 1e-3 times its mean eigenvalue (so that directions of
    % negligible energy, which make the canonical correlations of a nearly deterministic
    % series tie at 1, are not amplified)
    [evecs, evals] = eig(Yfut * Yfut');
    evals = max(diag(evals), 0);
    if ~(mean(evals) > 0)
        error('The series has no variation');
    end
    root = sqrt(evals + 1e-3 * mean(evals));
    W = (evecs ./ root') * evecs';
    [U, ~, ~] = svd(W * O);
    gamma = ((evecs .* root') * evecs') * U(:, 1:n); % extended observability matrix (r x n)

    gPinv = BF_Pinv(gamma(1:end-1, :));
    xNow = BF_Pinv(gamma) * O; % n x j
    xNext = gPinv * Ominus;
    C = gamma(1, :);
    A = gPinv * gamma(2:end, :);

    % Noise covariances from the residuals of the joint regression on the state
    z = [xNext; Yfut(1, :)];
    theta = z * BF_Pinv(xNow);
    rho = z - theta * xNow;
    Q = rho(1:n, :) * rho(1:n, :)' / j;
    R = rho(n+1, :) * rho(n+1, :)' / j;
    S = rho(1:n, :) * rho(n+1, :)' / j;

    A = SUB_stabilize(A);
    P = SUB_dare(A, C, Q, R, S, n);
    innov = C * P * C' + R;
    if any(~isfinite(P(:))) || ~(innov > 0)
        error('The Riccati equation has no usable solution');
    end
    K = (A * P * C' + S) / innov;

    % State signs: make every entry of C non-negative
    sgn = ones(n, 1); sgn(C < 0) = -1;
    A = A .* (sgn * sgn');
    K = K .* sgn;
    C = C .* sgn';

    yp = BF_SSPredict(A, K, C, y, 1);
    err = yp - y;
    err = err(~isnan(err));
    if isempty(err)
        fit.loss = NaN;
    else
        fit.loss = mean(err.^2);
    end
    fit.A = A; fit.K = K; fit.C = C; fit.order = n;
end
% ------------------------------------------------------------------------------
function A = SUB_stabilize(A)
    % Reflect eigenvalues with modulus above 1 to 1/conj(lambda)
    [V, D] = eig(A);
    lam = diag(D);
    unstable = abs(lam) > 1;
    if ~any(unstable)
        return
    end
    lam(unstable) = 1 ./ conj(lam(unstable));
    A = real(V * diag(lam) / V);
end
% ------------------------------------------------------------------------------
function P = SUB_dare(A, C, Q, R, S, n)
    % Stabilizing solution of the filter Riccati equation
    %   P = A P A' + Q - (A P C' + S)(C P C' + R)^(-1)(A P C' + S)'
    % by the structured doubling algorithm, after removing the cross term S.
    if ~(R > 0)
        error('The Riccati equation has no solution (R <= 0)');
    end
    Ak = A' - C' * (S' / R); % (the dual, control-form problem)
    Gk = C' * (C / R);
    Hk = Q - S * (S' / R);
    Hk = 0.5 * (Hk + Hk');
    I = eye(n);
    converged = false;
    for k = 1:100
        W = I + Gk * Hk;
        AW = Ak / W; % Ak * inv(W)
        Hnew = Hk + Ak' * Hk * (W \ Ak);
        Gk = Gk + AW * Gk * Ak';
        Ak = AW * Ak;
        Hnew = 0.5 * (Hnew + Hnew');
        Gk = 0.5 * (Gk + Gk');
        if any(~isfinite(Hnew(:)))
            break
        end
        delta = norm(Hnew - Hk, 1);
        Hk = Hnew;
        if delta <= 1e-13 * max(norm(Hk, 1), realmin)
            converged = true;
            break
        end
    end
    if ~converged
        error('The Riccati equation has no stabilizing solution');
    end
    P = Hk;
end
