function [A, K, C] = BF_ARMAtoSS(a, c)
% BF_ARMAtoSS   Innovations-form state-space model of an ARMA model (observer canonical form).
%
% Converts A(q) y(t) = C(q) e(t), with a = [1, a_1, ..., a_p] and c = [1, c_1, ..., c_q]
% the coefficients of A(q) = 1 + a_1 q^-1 + ... and C(q) = 1 + c_1 q^-1 + ..., to the
% model x(t+1) = A x(t) + K e(t), y(t) = C x(t) + e(t) of order n = max(p, q), in observer
% canonical form: the first column of A is -[a_1, ..., a_n], A has ones on its
% superdiagonal, K = [c_1 - a_1, ..., c_n - a_n] and C = [1, 0, ..., 0].
%
%---NOTES:
% L. Ljung, System Identification: Theory for the User, 2nd ed., Prentice Hall (1999),
% ch. 4; T. Kailath, Linear Systems, Prentice Hall (1980), ch. 2.

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

a = a(:)'; c = c(:)';
n = max(length(a), length(c)) - 1;
if n <= 0
    A = zeros(0); K = zeros(0, 1); C = zeros(1, 0);
    return
end
ap = zeros(1, n + 1); cp = zeros(1, n + 1);
ap(1:length(a)) = a;
cp(1:length(c)) = c;
A = zeros(n);
A(:, 1) = -ap(2:end)';
A(1:n-1, 2:n) = A(1:n-1, 2:n) + eye(n - 1);
K = (cp(2:end) - ap(2:end))';
C = [1, zeros(1, n - 1)];

end
