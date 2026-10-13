function [a, fpe, loss] = BF_ARfb(y, p)
% BF_ARfb   AR model fitted by forward-backward least squares.
%
% Fits A(q) y(t) = e(t), A(q) = 1 + a_1 q^-1 + ... + a_p q^-p, by minimizing the sum of
% the squared forward and backward (time-reversed) prediction errors jointly, using only
% the measured samples (no windowing, no mean removal): the modified covariance method.
%
%---INPUTS:
% y, the time series (length N > p)
% p, the model order (>= 1)
%
%---OUTPUTS:
% a, the polynomial coefficients [1, a_1, ..., a_p] (a_k is the negative of the usual AR
%       coefficient on lag k)
% fpe, Akaike's final prediction error, loss*(1 + p/N)/(1 - p/N)
% loss, the mean squared forward residual over the N - p samples that have p past values
%
%---NOTES:
% S. L. Marple Jr., Digital Spectral Analysis with Applications, Prentice-Hall (1987),
% ch. 8. FPE: H. Akaike, "Fitting autoregressive models for prediction", Ann. Inst.
% Statist. Math. 21, 243 (1969). The least-squares solution is the minimum-norm one
% (BF_Pinv, relative cutoff 1e-10) if the design is (nearly) rank deficient.

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
if N <= p
    error('The series is too short for an AR(%u) model', p);
end
H = hankel(y(1:p+1), y(p+1:N))'; % row j = y(j:j+p), j = 1..N-p
Xf = H(:, p:-1:1); tf = H(:, p+1); % forward: y(t) on y(t-1), ..., y(t-p)
Xb = H(:, 2:p+1); tb = H(:, 1); % backward: y(t) on y(t+1), ..., y(t+p)
D = [Xf; Xb];
phi = BF_Pinv(D) * [tf; tb];
res = tf - Xf * phi;
loss = (res' * res) / (N - p);
fpe = loss * (1 + p/N) / (1 - p/N);
a = [1, -phi'];

end
