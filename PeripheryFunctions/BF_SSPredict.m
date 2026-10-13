function yp = BF_SSPredict(A, K, C, y, steps)
% BF_SSPredict   k-step-ahead predictions of an innovations-form state-space model from a diffuse initial state.
%
% For the model x(t+1) = A x(t) + K e(t), y(t) = C x(t) + e(t) (unit innovation variance),
% the initial state is taken to be completely unknown (zero mean, infinite variance) and
% the exact diffuse Kalman filter is run: the covariance of the one-step state prediction
% is P* + kappa*Pinf with kappa -> infinity, Pinf starts as the identity, and the gain and
% update equations are the kappa -> infinity limits of the standard ones. The prediction
% yhat(t | t-k) = C A^(k-1) xhat(t-k+1) is returned only where its variance is finite; the
% other samples (the first n + k - 1 of an observable model of order n) are NaN.
%
%---INPUTS:
% A, the n x n state matrix
% K, the n-vector noise gain
% C, the n-vector output matrix
% y, the observed time series (length N)
% steps, the prediction horizon k >= 1 (default: 1)
%
%---OUTPUT:
% yp, the N x 1 vector of predictions (NaN where undefined)
%
%---NOTES:
% J. Durbin and S. J. Koopman, Time Series Analysis by State Space Methods, 2nd ed.,
% Oxford University Press (2012), ch. 5 (exact initial Kalman filtering), with the
% formulae adapted to the innovations form, in which the state disturbance K e(t) is
% correlated with the measurement error e(t). A model of order 0 predicts 0.

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

if nargin < 5 || isempty(steps)
    steps = 1;
end
y = y(:);
N = length(y);
n = size(A, 1);
if n == 0
    yp = zeros(N, 1);
    return
end
C = reshape(C, 1, n);
K = reshape(K, n, 1);
relTol = 1e-9;
row = C * A^(steps - 1); % maps the one-step predictor state to the k-step prediction

states = zeros(N, n);
dvar = zeros(N, 1); % diffuse part of the variance of row*x(t)
x = zeros(n, 1);
Pinf = eye(n);
Pst = zeros(n);
KK = K * K';
diffuse = true;
tolF = relTol * (C * C');
tolP = relTol * n;
for t = 1:N
    states(t, :) = x';
    v = y(t) - C * x;
    if diffuse
        dvar(t) = row * Pinf * row';
        Minf = Pinf * C';
        Mst = Pst * C';
        Finf = C * Minf;
        Fst = C * Mst + 1;
        Ust = A * Mst + K;
        if Finf > tolF
            Uinf = A * Minf;
            x = A * x + Uinf * (v / Finf);
            Pn = A * Pst * A' + KK - (Ust * Uinf' + Uinf * Ust') / Finf + (Uinf * Uinf') * (Fst / Finf^2);
            Pinf = A * Pinf * A' - (Uinf * Uinf') / Finf;
        else
            x = A * x + Ust * (v / Fst);
            Pn = A * Pst * A' + KK - (Ust * Ust') / Fst;
            Pinf = A * Pinf * A';
        end
        Pst = 0.5 * (Pn + Pn');
        Pinf = 0.5 * (Pinf + Pinf');
        if trace(Pinf) <= tolP
            diffuse = false;
        end
    else
        M = Pst * C';
        F = C * M + 1;
        U = A * M + K;
        x = A * x + U * (v / F);
        Pn = A * Pst * A' + KK - (U * U') / F;
        Pst = 0.5 * (Pn + Pn');
    end
end

yp = NaN(N, 1);
if N >= steps
    src = 1:(N - steps + 1);
    ok = dvar(src) <= relTol * max(row * row', realmin);
    vals = states(src, :) * row';
    vals(~ok) = NaN;
    yp(steps:end) = vals;
end

end
