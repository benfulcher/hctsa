function C = BF_NaNCov(X,makeCoeff,makeDist)
% BF_NaNCov     Covariance estimate including NaNs for an input matrix, X.
%
% Not exact, because removes full mean across all values, rather than across
% overlapping range, but should a reasonable approximation when number of NaNs
% is small.
%
% Output can be either the covariance matrix, or matrix of correlation
% coefficients, specified by the second input, makeCoeff.

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
% Check Inputs:
% ------------------------------------------------------------------------------

if nargin < 2
    makeCoeff = 1; % by default, convert to correlation coefficients
end

if nargin < 3
    makeDist = 0; % by default, don't convert to distances
end

% ------------------------------------------------------------------------------

if any(isnan(X(:)))
    % Indicate non-NaN values:
    goodValues = ~isnan(X);

    % Compute column means and standard deviations, excluding NaNs (vectorized:
    % exact match to a per-column loop, but without the per-column function-call
    % overhead, which matters once numCol reaches into the thousands, e.g.,
    % feature-feature correlations over many operations):
    colMeans = mean(X,1,'omitnan');
    colStds = std(X,0,1,'omitnan');

    % Scale each column to unit standard deviation (leaving columns with zero or
    % undefined spread unscaled), after removing its mean. This keeps the values
    % within single-precision range, whatever the scale of X: features above ~3.4e38
    % (and products of those above ~1.8e19) overflow single precision, which made
    % their correlations with every other feature infinite or undefined.
    colScale = colStds;
    isBadCol = ~(colScale > 0 & isfinite(colScale));
    colScale(isBadCol) = 1;
    Xc = (X - colMeans)./colScale;

    % X0 copies Xc but puts zeros over NaNs. Cast to single before the
    % O(numCol^2 * numRows) matrix multiply below, which dominates runtime
    % for large numCol: single precision is roughly 1.5-2x faster via BLAS (and
    % halves the memory of the result), at a relative error of ~1e-6, negligible
    % given this is already an approximate covariance estimate (see function
    % description above).
    X0 = single(Xc);
    X0(~goodValues) = 0; % NaN -> 0

    % Count good points (overlapping non-NaN values)
    goodValues = single(goodValues);
    goodBoth = goodValues' * goodValues;

    % This is our approximation to the covariance matrix of the scaled data, which is
    % the matrix of correlation coefficients:
    C = (X0' * X0) ./ (goodBoth - 1);

    if makeCoeff
        % Correlation with a column of zero or undefined spread is undefined:
        C(isBadCol,:) = NaN;
        C(:,isBadCol) = NaN;
    else
        % Undo the scaling to recover covariances (in double precision, as they can be
        % beyond the range of single precision):
        C = double(C).*(colScale'*colScale);
    end
else
    % no NaNs, use the matlab cov function:
    C = cov(X);

    if makeCoeff
        % Normalize by sample standard deviations:
        colStds = std(X,0,1);
        S = colStds'*colStds;
        C = C./S;
    end
end

% ------------------------------------------------------------------------------
% Convert to distances
% ------------------------------------------------------------------------------

if makeDist && makeCoeff
    C(logical(eye(size(C)))) = 1;
    C = 1 - C;
end

end
