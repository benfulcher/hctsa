function X = BF_Pinv(M, relTol)
% BF_Pinv   Pseudoinverse with a relative singular-value cutoff.
%
% The Moore-Penrose pseudoinverse of M from its singular value decomposition, keeping
% the singular values above relTol times the largest (default relTol = 1e-10). The
% cutoff is far above the level of rounding error, so that the numerical rank of a
% (nearly) rank-deficient matrix, such as the lagged design of an exactly predictable
% series, does not depend on rounding; singular values this small carry no statistical
% information for the least-squares fits that use this function.
%
%---INPUTS:
% M, a matrix
% relTol, the relative cutoff (default: 1e-10)
%
%---OUTPUT:
% X, the pseudoinverse (size(M') )

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

if nargin < 2 || isempty(relTol)
    relTol = 1e-10;
end
if isempty(M)
    X = zeros(size(M'));
    return
end
[U, S, V] = svd(M, 'econ');
sv = diag(S);
keep = sv > relTol * max(sv);
X = V(:, keep) * diag(1 ./ sv(keep)) * U(:, keep)';

end
