function R = BF_pdist(dataMatrix,distMetric,toVector,opts,beSilent,minPropGood)
% BF_pdist  Pairwise distances between rows of a data matrix.
%
% Same as pdist but then goes through and fills in NaNs with indiviually
% calculated values using an overlapping range of good values.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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
    toVector = 0;
end
if nargin < 4
    opts = [];
end
if nargin < 5
    beSilent = 0;
end
if nargin < 6
    % By default, don't require a minimum proportion of good values to be
    % present to compute a pairwise distance
    minPropGood = 0;
end

[n1, n2] = size(dataMatrix); % We're computing for rows (operations are rows)

% ------------------------------------------------------------------------------
% Define the distance function
% ------------------------------------------------------------------------------
switch distMetric
    case {'Euclidean','euclidean'}
        dij = @(v1,v2) sqrt(sum((v1-v2).^2))/length(v1)*n2; % if less entries, don't bias
    case {'corr','correlation','abscorr','Pearson'}
        dij = @(v1,v2) subdc(v1,v2,'Pearson');
    case 'Spearman'
        dij = @(v1,v2) subdc(v1,v2,'Spearman');
end

% ------------------------------------------------------------------------------
% Compute pairwise distances
% ------------------------------------------------------------------------------
switch distMetric
case 'mi'
    % Mutual information distances: can't make use of the inbuilt pdist function
    if ~isempty(opts)
        nbins = opts; % for MI, extra argument specifies nbins
    else
        nbins = 10;
    end
    if ~beSilent, fprintf(1,'Using a histogram with %u bins\n',nbins); end

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
            % mis(i,j) = BF_MutualInformation(dataMatrix(i,goodboth),dataMatrix(j,goodboth),'quantile','quantile',nbins); % by quantile with nbins
            mis(j,i) = mis(i,j);
        end
        if (mod(i,floor(n1/50)) == 0)
            fprintf(1,'Approximately %s remaining! We''re at %u / %u\n', ...
                        BF_thetime(toc(mitimer)/i*(n1-i)),i,n1);
        end
    end
    clear mitimer % stop timing
    R = mis; clear mis; % not really an R but ok.


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
    if ~beSilent, fprintf(1,' Done in %s.\n',BF_thetime(toc)); end


case {'euclidean','Euclidean','corr','correlation','abscorr'}
    % First use in-built pdist, which is fast
    if ~beSilent
        fprintf(1,'First computing pairwise distances using pdist...');
    end
    tic
    if strcmp(distMetric,'abscorr')
        R = pdist(dataMatrix,'corr');
    else
        R = pdist(dataMatrix,distMetric);
    end
    R = squareform(R); % Make a matrix
    if ~beSilent
        fprintf(1,' Done in %s.\n',BF_thetime(toc));
    end

    % Now go through and fill in any NaNs
    [nani, nanj] = find(isnan(R));
    if ~isempty(nani) % there are NaNs in R
        ij = (nanj >= nani); % only keep diagonal or upper diagonal entries
        nani = nani(ij);
        nanj = nanj(ij);
        NotNaN = ~isnan(dataMatrix);

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
                        BF_thetime(toc(NaNtimer)/i*(length(nani)-i)),i,length(nani));
            end
        end
        clear NaNtimer % stop the timer
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
    R = 1 - abs(1-R);
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

% ------------------------------------------------------------------------------
function d = subdc(v1,v2,corrType)
    rc = corr(v1,v2,'type',corrType);
    if isnan(rc)
        d = 2; % return maximum distance--for this subset, constant values
    else
        d = 1 - rc; % Return the (raw) correlation coefficient
    end
end
% ------------------------------------------------------------------------------

end
