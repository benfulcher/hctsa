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
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

% Number of rows and columns (should be the same):
[~,numCol] = size(X);

if any(isnan(X(:)))
    % Indicate non-NaN values:
    GoodValues = single(~isnan(X));

    % Compute column means, excluding NaNs:
    % Problem is with X(~isnan) -> 0
    meanNotNan = @(x) mean(x(~isnan(x)));
    ColMeans = arrayfun(@(x)meanNotNan(X(:,x)),1:numCol);

    % Remove mean from each column, to make centered version:
    Xc = bsxfun(@minus,X,ColMeans);

    % X0 copies Xc but puts zeros over NaNs:
    X0 = Xc;
    X0(~GoodValues) = 0; % NaN -> 0

    % Count good points (overlapping non-NaN values)
    GoodBoth = GoodValues' * GoodValues;

    % This is our approximation to the covariance matrix:
    C = (X0' * X0) ./ (GoodBoth - 1);

    % Convert to a correlation coefficient:
    if makeCoeff
        % Normalize by sample standard deviations:
        stdNotNan = @(x) std(x(~isnan(x)));
        ColStds = arrayfun(@(x)stdNotNan(X(:,x)),1:numCol);
        S = ColStds'*ColStds;
        C = C./S;
    end
else
    % no NaNs, use the matlab cov function:
    C = cov(X);

    if makeCoeff
        % Normalize by sample standard deviations:
        ColStds = arrayfun(@(x)std(X(:,x)),1:numCol);
        S = ColStds'*ColStds;
        C = C./S;
    end
end

% ------------------------------------------------------------------------------
% Convert to distances
% ------------------------------------------------------------------------------

if makeDist && makeCoeff
    C(logical(eye(size(C)))) = 1;
    C = 1-C;
end

end
