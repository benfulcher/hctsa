function R = BF_pdist(dataMatrix,distMetric,toVector,~,beSilent,minPropGood)
% BF_pdist  Pairwise distances between rows of a data matrix.
%
% Same as pdist but then goes through and fills in NaNs with values
% calculated using only the overlapping range of good values between each
% NaN-affected pair. For the correlation-based metrics (Pearson/Spearman/
% corr/correlation/abscorr) this NaN fix-up is a single vectorized
% corr(...,'Rows','pairwise') call covering every affected row at once
% (which reproduces a naive per-pair loop exactly); Euclidean distance
% has no native pairwise-complete equivalent, so it still uses an explicit
% per-pair loop.
%
% A correlation that is undefined for a pair (e.g., a constant row) is given the
% maximum distance for the metric: 2 for the signed correlation distances, 1 for
% 'abscorr'.
%
%---INPUTS:
% dataMatrix, the data matrix; distances are computed between its rows
% distMetric, the distance metric (default: 'euclidean'):
%       'euclidean'/'Euclidean': Euclidean distance (for pairs with missing values,
%                   computed on the overlapping good values and rescaled by
%                   sqrt(numFeatures/numGoodValues))
%       'Pearson'/'PearsonCorr'/'corr'/'correlation': 1 - Pearson correlation
%       'Spearman': 1 - Spearman correlation
%       'abscorr': 1 - |Spearman correlation|
%       'corr_fast'/'abscorr_fast': 1 - R and 1 - |R| using BF_NaNCov (an approximation
%                   to the Pearson correlation when there are missing values)
%       'mi-gaussian': Gaussian mutual information (nats) between rows, a similarity
%                   (not a distance; its diagonal is Inf and toVector is not possible).
%                   'mi' is the deprecated name for this metric.
% toVector, return the distances as a vector, as pdist does (default: false, a matrix)
% (the fourth input is unused; pass [])
% beSilent, suppress progress messages (default: false)
% minPropGood, the minimum proportion of overlapping good values required to compute
%       a distance for a pair (default: 0)

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
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(distMetric)
    distMetric = 'euclidean';
    fprintf(1,'Using the Euclidean distance metric\n');
end
if nargin < 3 || isempty(toVector)
    toVector = false;
end
if nargin < 5
    beSilent = false;
end
if nargin < 6
    % By default, don't require a minimum proportion of good values to be
    % present to compute a pairwise distance
    minPropGood = 0;
end

% We're computing for rows (operations are rows):
[n1, n2] = size(dataMatrix);

% ------------------------------------------------------------------------------
% Compute pairwise distances
% ------------------------------------------------------------------------------
switch distMetric
case {'mi-gaussian','mi'}
    % Gaussian mutual information (a function of the correlation only): can't make
    % use of the inbuilt pdist function. ('mi' is the deprecated name for this metric.)
    if strcmp(distMetric,'mi')
        warning('hctsa:deprecated',['BF_pdist metric ''mi'' is deprecated: it is the Gaussian ' ...
                    'mutual information, now named ''mi-gaussian''.']);
    end
    % Gaussian mutual information is -0.5*log(1-r^2) for correlation r, using the
    % overlapping good values of each pair:
    r = corr(dataMatrix','Rows','pairwise');
    R = -0.5*log(1 - r.^2);
    R(1:n1+1:end) = Inf; % a variable has unbounded information about itself

case {'corr_fast','abscorr_fast'}
    % Try using fast approximation to correlation coefficients when data includes NaNs
    % This is an approximation in that it centers columns on their full mean rather than
    % the mean of overlapping good values, but it's a start, and a good approximation
    % for a small proportion of NaNs.
    % Ben Fulcher, 2014-06-26
    if ~beSilent
        fprintf(1,'Using BF_NaNCov to approximate correlations between %u objects...',size(dataMatrix,1));
    end
    tic
    R = BF_NaNCov(dataMatrix',1,1);
    if ~beSilent
        fprintf(1,' Done in %s.\n',BF_TheTime(toc));
    end

case {'euclidean','Euclidean','PearsonCorr','Pearson','Spearman','corr','correlation','abscorr'}

    % Distance for an undefined correlation (before the abscorr transform at the end,
    % which maps R = 1 - rho to 1 - |rho|, so rho = 0 is the largest distance there):
    if strcmp(distMetric,'abscorr')
        undefinedDist = 1;
    else
        undefinedDist = 2;
    end

    isCorrType = ~ismember(distMetric,{'Euclidean','euclidean'});

    if isCorrType && ~any(isnan(dataMatrix(:)))
        % Fast path: no NaNs anywhere, so there is nothing for the fixup
        % below to ever do -- skip pdist+squareform entirely (which carry
        % nontrivial fixed overhead of their own even with nothing to fix)
        % and compute the full distance matrix in one vectorized call.
        if ismember(distMetric,{'Spearman','abscorr'})
            corrType = 'Spearman'; % matches the 'Spearman' forced for abscorr below
        else
            corrType = 'Pearson'; % matches pdist's 'corr'/'correlation' convention
        end
        if ~beSilent
            fprintf(1,'No NaNs present -- computing %s correlation distances directly...',corrType);
        end
        tic
        R = 1 - corr(dataMatrix','Type',corrType);
        R(isnan(R)) = undefinedDist; % undefined (e.g., constant-value row) correlation -> max distance
        R(1:n1+1:end) = 0; % exact self-distance (matches squareform's convention; applied last, overriding the fallback above for a degenerate row's self-entry)
        if ~beSilent
            fprintf(1,' Done in %s.\n',BF_TheTime(toc));
        end
        nani = []; % nothing left to fix up below
    else
        % ------------------------------------------------------------------------------
        % Define the distance function (Euclidean only -- the correlation-based
        % metrics are fixed up via a vectorized corr(...,'Rows','pairwise') call
        % below instead, so no per-pair distance function is needed for them)
        % ------------------------------------------------------------------------------
        if ismember(distMetric,{'Euclidean','euclidean'})
            dij = @(v1,v2) sqrt(sum((v1-v2).^2)*n2/length(v1)); % if less entries, don't bias
        end

        % First use in-built pdist, which is fast
        if ~beSilent
            fprintf(1,'First computing pairwise distances using pdist...');
        end
        tic
        if strcmp(distMetric,'abscorr')
            % Feature-feature Correlations are Spearman by default
            R = pdist(dataMatrix,'Spearman');
        elseif ismember(distMetric,{'Pearson','PearsonCorr'})
            R = pdist(dataMatrix,'correlation');
        else
            R = pdist(dataMatrix,distMetric);
        end
        R = squareform(R); % Make a matrix
        if ~beSilent
            fprintf(1,' Done in %s.\n',BF_TheTime(toc));
        end

        % Now go through and fill in any NaNs
        [nani, nanj] = find(isnan(R));
    end
    if ~isempty(nani) % there are NaNs in R
        NotNaN = ~isnan(dataMatrix);
        % A row is "intrinsically" affected (missing data, or -- even with no
        % missing data -- zero variance/all-tied values giving an undefined
        % correlation) iff EVERY one of its off-diagonal entries in R is NaN:
        % pdist computes each row's own mean/rank once, so any such problem
        % propagates to that row's distance to every other row, not just rows
        % it happens to share missing positions with. This deliberately does
        % NOT use unique(nani) from find(isnan(R)) directly: since R is
        % symmetric, that would incorrectly flag every OTHER row too, merely
        % because its one NaN entry (against the genuinely bad row) shows up
        % as a "NaN row" reference for that innocent row as well.
        affectedRows = find(sum(isnan(R),2) >= n1-1);

        if ismember(distMetric,{'Euclidean','euclidean'})
            % No native pairwise-complete equivalent for this (overlap-rescaled)
            % Euclidean distance, so fall back to an explicit per-pair loop:
            ij = (nanj >= nani); % only keep diagonal or upper diagonal entries
            nani = nani(ij);
            nanj = nanj(ij);

            if ~beSilent
                fprintf(1,['Recalculating distances individually for %u NaN ' ...
                                'entries in the distance matrix...\n'],length(nani));
            end

            NaNtimer = tic; % time it
            for i = 1:length(nani)
                ii = nani(i);
                jj = nanj(i);
                goodboth = (NotNaN(ii,:) & NotNaN(jj,:));
                if mean(goodboth) > minPropGood
                    R(ii,jj) = dij(dataMatrix(ii,goodboth)',dataMatrix(jj,goodboth)'); % Calculate the distance
                else
                    R(ii,jj) = NaN; % Not enough good, overlapping set of values -- store as NaN.
                end
                R(jj,ii) = R(ii,jj); % Add the symmetrized entry

                % Give update on time remaining after 1000 iterations (if more than 10000 total iterations)
                % and then 5 more times...
                if ~beSilent && ((i==1000 && length(nani) > 10000) || (mod(i,floor(length(nani)/5))==0))
                    fprintf(1,'Approximately %s remaining! We''re at %u / %u.\n', ...
                            BF_TheTime(toc(NaNtimer)/i*(length(nani)-i)),i,length(nani));
                end
            end
            clear NaNtimer % stop the timer
        else
            % Correlation-based metrics: fix every NaN-affected row in a single
            % vectorized pairwise-complete correlation call (corr(...,'Rows',
            % 'pairwise') reproduces the per-pair loop's output exactly, since
            % both use the same overlapping-good-values convention) rather than
            % looping pair-by-pair -- while still only touching the rows that
            % actually need it, keeping the fast path fast when few rows are
            % NaN-affected.
            if ismember(distMetric,{'Spearman','abscorr'})
                corrType = 'Spearman';
            else
                corrType = 'Pearson'; % as for pdist's 'corr'/'correlation'
            end

            if ~beSilent
                fprintf(1,['Recalculating distances for %u NaN-affected row(s) using ' ...
                                'pairwise-complete %s correlation...'],length(affectedRows),corrType);
            end
            NaNtimer = tic;

            Csub = corr(dataMatrix(affectedRows,:)',dataMatrix','Type',corrType,'Rows','pairwise');
            Dsub = 1 - Csub;
            Dsub(isnan(Csub)) = undefinedDist; % undefined (e.g., constant-value) correlation -> max distance
            overlapProp = (NotNaN(affectedRows,:) * NotNaN') / n2; % proportion of overlapping good values, per pair
            Dsub(overlapProp <= minPropGood) = NaN; % not enough good, overlapping set of values

            R(affectedRows,:) = Dsub;
            R(:,affectedRows) = Dsub';
            % squareform's diagonal convention is always exactly 0 (self-
            % distance is never computed by pdist in the first place); restore
            % that here since the vectorized self-correlation above is
            % otherwise undefined (NaN, mapped to 2 above) for a degenerate
            % (e.g. constant-valued) affected row correlated with itself:
            R(1:n1+1:end) = 0;

            if ~beSilent
                fprintf(1,' Done in %s.\n',BF_TheTime(toc(NaNtimer)));
            end
            clear NaNtimer
        end

        if any(isnan(R(:)))
            warning('%u pairs still produce NaNs, with less than %.3f%% overlap',sum(isnan(R(:)))/2,minPropGood);
        end
    end
otherwise
    error('Unknown distance metric ''%s''',distMetric);
end

% ------------------------------------------------------------------------------
% Transform from correlation distance to absolute correlation distance:
% ------------------------------------------------------------------------------
if ismember(distMetric,{'abscorr','abscorr_fast'})
    R = 1 - abs(1 - R);
end

% ------------------------------------------------------------------------------
% Convert from matrix back to a vector:
% ------------------------------------------------------------------------------
if toVector
    try
        R = squareform(R); % back to vector
    catch
        error('This metric is not consistant with a pairwise distance matrix...?')
    end
end

end
