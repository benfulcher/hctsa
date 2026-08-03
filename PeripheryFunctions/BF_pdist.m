function R = BF_pdist(dataMatrix,distMetric,toVector,opts,beSilent,minPropGood)
% BF_pdist  Pairwise distances between rows of a data matrix.
%
% Same as pdist but then goes through and fills in NaNs with values
% calculated using only the overlapping range of good values between each
% NaN-affected pair. For the correlation-based metrics (Pearson/Spearman/
% corr/correlation/abscorr) this NaN fix-up is a single vectorized
% corr(...,'Rows','pairwise') call covering every affected row at once
% (verified to reproduce a naive per-pair loop exactly); Euclidean distance
% has no native pairwise-complete equivalent, so it still uses an explicit
% per-pair loop.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
if nargin < 4
    opts = [];
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
case 'mi'
    % Mutual information distances: can't make use of the inbuilt pdist function
    if ~isempty(opts)
        numBins = opts; % for MI, extra argument specifies numBins
    else
        numBins = 10;
    end
    if ~beSilent,
        fprintf(1,'Using a histogram with %u bins\n',numBins);
    end

    goodies = ~isnan(dataMatrix); % now we can deal with NaNs into design matrix

    mis = zeros(n1);
    mitimer = tic; % Way faster to not store the time taken for every iteration
    for i = 1:n1
        % tic
        goodi = goodies(i,:);
        for j = i:n1
            goodj = goodies(j,:);
            goodboth = (goodi & goodj);
            % Using Information Dynamics Toolkit:
            mis(i,j) = IN_MutualInfo(dataMatrix(i,goodboth),dataMatrix(j,goodboth),'gaussian');
            % mis(i,j) = BF_MutualInformation(dataMatrix(i,goodboth),dataMatrix(j,goodboth),'quantile','quantile',numBins); % by quantile with numBins
            mis(j,i) = mis(i,j);
        end
        if (mod(i,floor(n1/50)) == 0)
            fprintf(1,'Approximately %s remaining! We''re at %u / %u\n', ...
                        BF_TheTime(toc(mitimer)/i*(n1-i)),i,n1);
        end
    end
    R = mis; % not really an R but ok.
    clear('mitimer','mis');

case {'corr_fast','abscorr_fast'}
    % Try using fast approximation to correlation coefficients when data includes NaNs
    % This is an approximation in that it centers columns on their full mean rather than
    % the mean of overlapping good values, but it's a start, and a good approximation
    % for a small proportion of NaNs.
    % Ben Fulcher, 2014-06-26
    if ~beSilent,
        fprintf(1,'Using BF_NaNCov to approximate correlations between %u objects...',size(dataMatrix,1));
    end
    tic
    R = BF_NaNCov(dataMatrix',1,1);
    if ~beSilent
        fprintf(1,' Done in %s.\n',BF_TheTime(toc));
    end

case {'euclidean','Euclidean','PearsonCorr','Pearson','Spearman','corr','correlation','abscorr'}

    isCorrType = ~ismember(distMetric,{'Euclidean','euclidean','PearsonCorr','Pearson'});

    if isCorrType && ~any(isnan(dataMatrix(:)))
        % Fast path: no NaNs anywhere, so there is nothing for the fixup
        % below to ever do -- skip pdist+squareform entirely (which carry
        % nontrivial fixed overhead of their own even with nothing to fix:
        % measured ~2x slower than going straight to corr() on a 500x1000
        % matrix) and compute the full distance matrix in one vectorized call.
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
        R(isnan(R)) = 2; % subdc's fallback: undefined (e.g. constant-value row) correlation -> max distance
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
            dij = @(v1,v2) sqrt(sum((v1-v2).^2))/length(v1)*n2; % if less entries, don't bias
        end

        % First use in-built pdist, which is fast
        if ~beSilent
            fprintf(1,'First computing pairwise distances using pdist...');
        end
        tic
        if strcmp(distMetric,'abscorr')
            % Feature-feature Correlations are Spearman by default
            R = pdist(dataMatrix,'Spearman');
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
            % Euclidean distance, so fall back to the original per-pair loop:
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
            % NaN-affected (see PeripheryFunctions/BF_pdist.m investigation
            % notes/commit message for benchmarks across NaN-density regimes).
            if ismember(distMetric,{'PearsonCorr','Pearson'})
                corrType = 'Pearson';
            else
                corrType = 'Spearman';
            end

            if ~beSilent
                fprintf(1,['Recalculating distances for %u NaN-affected row(s) using ' ...
                                'pairwise-complete %s correlation...'],length(affectedRows),corrType);
            end
            NaNtimer = tic;

            Csub = corr(dataMatrix(affectedRows,:)',dataMatrix','Type',corrType,'Rows','pairwise');
            Dsub = 1 - Csub;
            Dsub(isnan(Csub)) = 2; % subdc's fallback: undefined (e.g. constant-value) correlation -> max distance
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
