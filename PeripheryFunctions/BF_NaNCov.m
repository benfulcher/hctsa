% Covariance including NaNs for an input matrix, X
% Not exact, because removes full mean across all values, rather than across
% overlapping range, but should a reasonable approximation when number of NaNs
% is small.
% Output can be either the covariance matrix, or matrix of correlation
% coefficients, depending on the second input.
% ------------------------------------------------------------------------------
% Ben Fulcher, 2014-06-26
% ------------------------------------------------------------------------------

function C = NaNCov(X,MakeCoeff,MakeDist)

if nargin < 2
    MakeCoeff = 1; % by default, convert to correlation coefficients
end

if nargin < 3
    MakeDist = 0; % by default, don't convert to distances
end

% Number of rows and columns (should be the same):
[NRow,NCol] = size(X);

if any(isnan(X(:)))
    % Indicate non-NaN values:
    GoodValues = single(~isnan(X));
    
    % Compute column means, excluding NaNs:
    % Problem is with X(~isnan) -> 0
    meanNotNan = @(x) mean(x(~isnan(x)));
    ColMeans = arrayfun(@(x)meanNotNan(X(:,x)),1:NCol);

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
    if MakeCoeff
        % Normalize by sample standard deviations:
        stdNotNan = @(x) std(x(~isnan(x)));
        ColStds = arrayfun(@(x)stdNotNan(X(:,x)),1:NCol);
        S = ColStds'*ColStds;
        C = C./S;
    end
else
    % no NaNs, use the matlab cov function:
    C = cov(X);
    
    if MakeCoeff
        % Normalize by sample standard deviations:
        ColStds = arrayfun(@(x)std(X(:,x)),1:NCol);
        S = ColStds'*ColStds;
        C = C./S;
    end
end

% ------------------------------------------------------------------------------
% Convert to distances
% ------------------------------------------------------------------------------

if MakeDist && MakeCoeff
    C(logical(eye(size(C)))) = 1;
    C = 1-C;
end

end