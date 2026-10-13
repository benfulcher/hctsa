function [a, c, coefSD, loss] = BF_ARMAfit(y, p, q)
% BF_ARMAfit   ARMA(p, q) model fitted by conditional least squares from a Hannan-Rissanen start.
%
% Fits A(q) y(t) = C(q) e(t), with A(q) = 1 + a_1 q^-1 + ... + a_p q^-p and
% C(q) = 1 + c_1 q^-1 + ... + c_q q^-q (zero mean), by minimizing the conditional sum of
% squares V = sum_{t > r} e(t)^2, r = max(p, q), where e(t) = 0 for t <= r and
% e(t) = y(t) + a_1 y(t-1) + ... + a_p y(t-p) - c_1 e(t-1) - ... - c_q e(t-q) for t > r.
%
% (1) Start: the Hannan-Rissanen two-stage regression. A long AR model of order
% m = max(p + q, min(20, floor(N/10))) is fitted by least squares to y(m+1:N) and its
% residuals estimate the innovations; y(t) is then regressed on y(t-1), ..., y(t-p) and
% the estimated innovations at lags 1, ..., q, for t = m + max(p, q) + 1, ..., N. Zeros of
% C on or outside the unit circle are reflected to 1/conj(z).
% (2) Refinement: Levenberg-Marquardt iterations on V (step (J'J + lambda*diag(J'J))
% \ (-J'e), with lambda starting at 1e-2, divided by 10 after a step that lowers V and
% multiplied by 10 after one that does not, or that would make C non-invertible; at most
% 100 steps; stopping when a step lowers V by no more than 1e-12*V or lambda exceeds
% 1e10).
%
%---INPUTS:
% y, the time series
% p, q, the AR and MA orders
%
%---OUTPUTS:
% a, [1, a_1, ..., a_p] (a_k is the negative of the usual AR coefficient on lag k)
% c, [1, c_1, ..., c_q] (all zeros inside the unit circle)
% coefSD, the standard errors of [a_1..a_p, c_1..c_q]: the square roots of the diagonal
%       of s2*pinv(J'J) at the optimum, with J the Jacobian of the residuals and
%       s2 = V/(N - r - p - q)
% loss, the mean squared one-step prediction error of the fitted model on y, from a
%       diffuse initial state (BF_SSPredict; the first max(p, q) samples are not predicted)
%
%---NOTES:
% E. J. Hannan and J. Rissanen, "Recursive estimation of mixed autoregressive-moving
% average order", Biometrika 69, 81 (1982); P. J. Brockwell and R. A. Davis, Introduction
% to Time Series and Forecasting, 3rd ed., Springer (2016), sec. 5.1.4 (Hannan-Rissanen)
% and sec. 5.2 (conditional least squares); D. W. Marquardt, "An algorithm for
% least-squares estimation of nonlinear parameters", J. Soc. Indust. Appl. Math. 11, 431
% (1963). The least-squares solutions of the start use BF_Pinv (relative singular-value
% cutoff 1e-10), as does the covariance.

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
if q > 0
    m = max(p + q, min(20, floor(N/10)));
else
    m = 0;
end
t0 = m + max(p, q); % last sample not used as a target in stage 2
n2 = N - t0;
if n2 <= p + q + 1 || N - m <= m + 1
    error('The series is too short for an ARMA(%u,%u) model', p, q);
end

% ------------------------------------------------------------------------------
%% (1) Hannan-Rissanen start
% ------------------------------------------------------------------------------
% Stage 1: innovations from a long AR model
e = NaN(N, 1);
if m > 0
    H = hankel(y(1:m+1), y(m+1:N))';
    X1 = H(:, m:-1:1);
    t1 = H(:, m+1);
    phiLong = BF_Pinv(X1) * t1;
    e(m+1:N) = t1 - X1 * phiLong;
end
% Stage 2: regression on lagged values and lagged innovations
tt = (t0 + 1:N)';
X2 = zeros(n2, p + q);
for k = 1:p
    X2(:, k) = y(tt - k);
end
for k = 1:q
    X2(:, p + k) = e(tt - k);
end
beta = BF_Pinv(X2) * y(tt);
a = [1, -beta(1:p)'];
c = SUB_reflect([1, beta(p+1:end)']);

% ------------------------------------------------------------------------------
%% (2) Levenberg-Marquardt refinement of the conditional sum of squares
% ------------------------------------------------------------------------------
r = max(p, q);
theta = [a(2:end), c(2:end)]';
[res, J] = SUB_residuals(theta);
V = res' * res;
lambda = 1e-2;
for iter = 1:100*(V > 0) % (no refinement of an exact fit)
    g = J' * res;
    JJ = J' * J;
    D = diag(max(diag(JJ), realmin));
    accepted = false;
    while lambda <= 1e10
        delta = -(JJ + lambda * D) \ g;
        if any(~isfinite(delta))
            lambda = lambda * 10; % (a singular system)
            continue
        end
        thetaNew = theta + delta;
        if q > 0 && any(abs(roots([1, thetaNew(p+1:end)'])) >= 1)
            lambda = lambda * 10; % keep C invertible
            continue
        end
        resNew = SUB_residuals(thetaNew);
        Vnew = resNew' * resNew;
        if Vnew < V
            accepted = true;
            break
        end
        lambda = lambda * 10;
    end
    if ~accepted
        break
    end
    decrease = V - Vnew;
    theta = thetaNew;
    [res, J] = SUB_residuals(theta);
    V = res' * res;
    lambda = max(lambda / 10, 1e-12);
    if decrease <= 1e-12 * V
        break
    end
end
a = [1, theta(1:p)'];
c = [1, theta(p+1:end)'];

% Standard errors: s2 * diag(pinv(J'J)), from the singular value decomposition of J (the
% cutoff of BF_Pinv), so that the condition number is not squared
s2 = V / (N - r - p - q);
[~, Sv, Vv] = svd(J, 'econ');
sv = diag(Sv);
keep = sv > 1e-10 * max(sv);
coefSD = sqrt(s2 * sum((Vv(:, keep) ./ sv(keep)').^2, 2))';

if nargout > 3
    [A, K, C] = BF_ARMAtoSS(a, c);
    err = BF_SSPredict(A, K, C, y, 1) - y;
    err = err(~isnan(err));
    loss = mean(err.^2);
end

% ------------------------------------------------------------------------------
function [res, Jac] = SUB_residuals(th)
    % Conditional residuals e(r+1:N) (and their Jacobian with respect to th)
    aa = [1, th(1:p)'];
    cc = [1, th(p+1:end)'];
    w = filter(aa, 1, y);
    w(1:r) = 0;
    ee = filter(1, cc, w);
    res = ee(r+1:N);
    if nargout > 1
        Jac = zeros(N - r, p + q);
        for kk = 1:p
            v = [zeros(kk, 1); y(1:N-kk)];
            v(1:r) = 0;
            z = filter(1, cc, v);
            Jac(:, kk) = z(r+1:N);
        end
        for kk = 1:q
            u = [zeros(kk, 1); ee(1:N-kk)];
            z = -filter(1, cc, u);
            Jac(:, p + kk) = z(r+1:N);
        end
    end
end
% ------------------------------------------------------------------------------
function cOut = SUB_reflect(cIn)
    % Reflect zeros of C on or outside the unit circle to 1/conj(z)
    cOut = cIn;
    if length(cIn) < 2
        return
    end
    z = roots(cIn);
    bad = abs(z) >= 1;
    if any(bad)
        z(bad) = 1 ./ conj(z(bad));
        cOut = real(poly(z));
    end
end

end
