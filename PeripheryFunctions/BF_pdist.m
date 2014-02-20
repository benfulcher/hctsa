% ------------------------------------------------------------------------------
% BF_pdist
% ------------------------------------------------------------------------------
% 
% Same as pdist but then goes through and fills in NaNs with indiviually
% calculated values using an overlapping range of good values.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function R = BF_pdist(F,DistMetric,ToVector,opts)

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(DistMetric)
    DistMetric = 'euclidean';
    fprintf(1,'Using the Euclidean distance metric\n')
end
if nargin < 3 || isempty(ToVector)
    ToVector = 0;
end
if nargin < 4
    opts = [];
end

% if ~ismember(DistMetric,{'euclidean','corr','correlation','abscorr'})
%     error(['I''ve only done Euclidean and correlation so far you know, not ' DistMetric '!! SORRY!!'])
% end

[n1, n2] = size(F); % We're computing for rows (operations are rows)

% Define the distance function
switch DistMetric
    case {'Euclidean','euclidean'}
        dij = @(v1,v2) sqrt(sum((v1-v2).^2))/length(v1)*n2; % if less entries, don't bias
    case {'corr','correlation','abscorr'}
        dij = @(v1,v2) subdc(v1,v2);
end

if strcmp(DistMetric,'mi')
    % Mutual information distances: can't make use of the inbuilt pdist function
    if ~isempty(opts)
        nbins = opts; % for MI, extra argument specifies nbins
    else
        nbins = 10;
    end
    fprintf(1,'Using a histogram with %u bins\n',nbins)

    goodies = ~isnan(F); % now we can deal with NaNs into design matrix

    mis = zeros(n1);
    % times = zeros(n1,1);
    timer = tic; % Way faster to not store the time taken for every iteration
    for i = 1:n1
        % tic
        goodi = goodies(i,:);
        for j = i:n1
            goodj = goodies(j,:);
            goodboth = (goodi & goodj);
            mis(i,j) = benmi(F(i,goodboth),F(j,goodboth),'quantile','quantile',nbins); % by quantile with nbins
            mis(j,i) = mis(i,j);
        end
        % times(i) = toc;
        if (mod(i,floor(n1/50)) == 0)
            fprintf(1,'Approximately %s remaining! We''re at %u / %u\n', ...
                        BF_thetime(toc(timer)/i*(n1-i)),i,n1)
%             save(fn,'mis','m_ids_keep','-v7.3');
%             disp('saved')
        end
    end
    clear timer % stop timing
    R = mis; clear mis; % not really an R but ok.
else
    % First use in-built pdist, which is fast
    fprintf(1,'First computing pairwise distances using pdist...');
    tic
    if strcmp(DistMetric,'abscorr')
        R = pdist(F,'corr');
    else
        R = pdist(F,DistMetric);
    end
    R = squareform(R); % Make a matrix
    fprintf(1,' Done in %s.\n',BF_thetime(toc));
    
    % Now go through and fill in any NaNs
    [nani, nanj] = find(isnan(R));
    if ~isempty(nani) % there are NaNs in R
        ij = (nanj>=nani); % only keep diagonal or upper diagonal entries
        nani = nani(ij);
        nanj = nanj(ij);
        goodies = ~isnan(F);
        
        fprintf(1,'Recalculating distances individually for %u NaN entries in the distance matrix...\n',length(nani));
        
        % Storing individual times for large numbers of calculations can be inefficient
        % times = zeros(length(nani),1); 
        timer = tic;
        for i = 1:length(nani)
            ii = nani(i);
            jj = nanj(i);
            goodi = goodies(ii,:);
            goodj = goodies(jj,:);
            goodboth = (goodi & goodj);
            R(ii,jj) = dij(F(ii,goodboth),F(jj,goodboth)); % calculate the distance
            R(jj,ii) = R(ii,jj); % add the symmetrized entry
            % times(i) = toc;
            
            % Give update on time remaining after 10 iterations (if more than 100 total iterations)
            % and then 10 more times...
            if (i==10 && length(nani) > 100) || (mod(i,floor(length(nani)/10))==0)
                fprintf(1,'Approximately %s remaining! We''re at %u / %u\n', ...
                        BF_thetime(toc(timer)/i*(length(nani)-i)),i,length(nani))
            end
        end
        clear timer % stop the timer
    end
end

if strcmp(DistMetric,'abscorr')
    % Transformation from correlation to absolute correlation distance:
    R = 1 - abs(1-R);
end


if ToVector
    try
        R = squareform(R); % back to vector
    catch
        error('This metric is not consistant with a pairwise distance matrix...?')
    end
end

    function d = subdc(v1,v2)
        rc = corrcoef(v1,v2);
        if isnan(rc), d = 2; % return maximum distance--for this subset, constant values
        else d = 1-rc(2,1); % return the (raw) correlation coefficient
        end
    end

end