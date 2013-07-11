function out = SD_surrogatetest(x,surrmeth,nsurrs,extrap,teststat)
% Ben Fulcher 28/1/2011
% Looks at test statistic of surrogates compared to that of the given time
% series
% first four inputs are inputs to SD_makesurrogates, teststatistic it the
% teststatistic to compare -- can specify many in a cell, will return
% outputs for each teststat specified. Better this way because creating
% surrogates is the expensive step.

%% CHECK INPUTS

if nargin < 2 || isempty(surrmeth)
    surrmeth = 'RP'; % Randomize phases
end
if nargin < 3 || isempty(nsurrs)
    nsurrs = 99; % create 99 surrogates for a 0.01 significance level 1-sided test
end
if nargin < 4 || isempty(nsurrs)
    extrap = [];
end
if nargin < 5 || isempty(teststat)
    teststat = 'AMI'; % automutual information
end

if ischar(teststat)
    teststat = {teststat};
end

N = length(x);

%% OK NOW DO SHIT

% 1) make surrogates
z = SD_makesurrogates(x,surrmeth,nsurrs,extrap);
% z is matrix where each column is a surrogate time series

% 2) evaluate test statistic on each surrogate
if ismember('ami1',teststat)
    % look at AMI(1) of surrogates compared to that of signal itself
    % This statistic is used by Nakamura et al. (2006), PRE
    % could use CO_ami_benhist or TSTL, but I'll use benmi
    % Apparently there are upper and lower bounds on the number of bins to
    % use: [1+log_2(N)], [sqrt(N)]
    nbins = ceil(1+log2(N)); %round(mean([1+log2(N),sqrt(N)]));
    AMIx = benmi(x(1:end-1),x(2:end),'quantile','quantile',nbins);
    AMIsurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        AMIsurr(i) = benmi(z(1:end-1,i),z(2:end,i),'quantile','quantile',nbins);
    end
    % so we have a value AMIx, and a distribution for the surrogates
    % AMIsurr -- we must compare and return something meaningful
    % surrogates should always have lower AMI than original signal
    somestats = SDgivemestats(AMIx,AMIsurr,'right');
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.ami_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end

if ismember('fmmi',teststat)
    % look at first minimum of mutual information of surrogates compared to
    % that of signal itself
    fmmix = CO_fmmi(x);
    fmmisurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        try
            fmmisurr(i) = CO_fmmi(z(i,:));
        catch
            out = NaN; return
        end
    end
    if any(isnan(fmmisurr))
        RA_keyboard
    end
    
    % FMMI should be higher for signal than surrogates
    somestats = SDgivemestats(fmmix,fmmisurr,'right');
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.fmmi_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end


if ismember('o3',teststat) % third order statistic in Schreiber, Schmitz (Physica D)
    tau = 1;
    o3x = 1/(N-tau)*sum((x(1+tau:end) - x(1:end-tau)).^3);
    o3surr = zeros(nsurrs,1);
    for i = 1:nsurrs
        o3surr(i) = 1/(N-tau)*sum((z(1+tau:end,i) - z(1:end-tau,i)).^3);
    end
    somestats = SDgivemestats(o3x,o3surr,'both');
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.o3_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end

if ismember('tc3',teststat) % TC3 statistic -- another time-reversal asymmetry measure
    tau = 1;
    tmp = CO_tc3(x,tau);
    tc3x = tmp.raw;
    tc3surr = zeros(nsurrs,1);
    for i = 1:nsurrs
        tmp = CO_tc3(z(:,i),tau);
        tc3surr(i) = tmp.raw;
    end
    
    somestats = SDgivemestats(tc3x,tc3surr,'both');
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.tc3_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end

if ismember('nlpe',teststat) % locally constant phase space prediction error
    fprintf(1,'''nlpe'' can be very time consuming\n')
    de = 3; tau = 1; % embedding parameters: fixed like a dummy!
    tmp = MS_nlpe(x,de,tau);
    nlpex = tmp.msqerr;
    nlpesurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        tmp = MS_nlpe(z(:,i),de,tau);
        nlpesurr(i) = tmp.msqerr;
    end
    
    somestats = SDgivemestats(nlpex,nlpesurr,'right'); % NLPE should be higher than surrogates
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.nlpe_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end

if ismember('fnn',teststat)
    fprintf(1,'fnn takes forever...\n')
    % false nearest neighbours at d=2;
    tmp = MS_fnn(x,2,1,5,1);
    fnnx = tmp.pfnn_2;
    fnnsurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        tmp = MS_fnn(z(:,i),2,1,5,1);
        fnnsurr(i) = tmp.pfnn_2;
    end
    
    somestats = SDgivemestats(fnnx,fnnsurr,'right'); % FNN(2) should be higher than surrogates?
    fnames = fieldnames(somestats);
    for i = 1:length(fnames)
        eval(sprintf('out.fnn_%s = somestats.%s;',fnames{i},fnames{i}));
    end
end





function somestats = SDgivemestats(statx,statsurr,leftrightboth)
    nsurrs = length(statsurr);
    if any(isnan(statsurr))
        RA_keyboard
    end
%         leftrightboth = {'left','right','both'}
    % have a distribution of some statistic with samples statsurr
    % and we have the value of the statistic for a given process in
    % statx. Want to return measures of how consistant the measured
    % statistic is with the sample statsurr.

    % ASSUME GAUSSIAN DISTRIBUTION:
    % so can use 1/2-sided z-statistic
    [~,p,~,zscore] = ztest(statx,mean(statsurr),std(statsurr),0.05,leftrightboth);
    somestats.p = p; % pvalue
    somestats.zscore = zscore; % z-statistic

    % fit a kernel distribution to zscored distribution:
    if std(statsurr) == 0
        % all surrogates have same value of this statistic
        % cannot do a meaningful zscore -- do it raw
        [f, xi] = ksdensity(statsurr);
        % find where the statx value would be:
        if statx < min(xi) || statx > max(xi)
            somestats.f = 0; % out of range -- assume p=0 here
        else
            [~, minhere] = min(abs(statx-xi));
            somestats.f = f(minhere); % return probability density where the point is
        end    
    else
        zscstatsurr = (statsurr-mean(statsurr))/std(statsurr);
        zscstatx = (statx-mean(statsurr))/std(statsurr);
        [f, xi] = ksdensity(zscstatsurr);

        % find where the statx value would be:
        if zscstatx < min(xi) || zscstatx > max(xi)
            somestats.f = 0; % out of range -- assume p=0 here
        else
            [~, minhere] = min(abs(zscstatx-xi));
            somestats.f = f(minhere); % return probability density where the point is
        end
    end
    
    % what fraction of the range is the sample in?
    medsurr = median(statsurr);
    iqrsurr = iqr(statsurr);
    if iqrsurr == 0
        somestats.mediqr = NaN;
    else
        somestats.mediqr = abs(statx-medsurr)/iqrsurr;
    end

    % rank statistic
    [~, ix] = sort([statx; statsurr]);
    xfitshere = find(ix == 1)-1;
    if strcmp(leftrightboth,'right') % x statistic smaller than distribution
        xfitshere = nsurrs + 1 - xfitshere; % how far from highest ranked value
    elseif strcmp(leftrightboth,'both')
        xfitshere = min(xfitshere,nsurrs+1-xfitshere);
    end

    if isempty(xfitshere)
        somestats.prank = 1/(nsurrs+1); % rank-based p-value
    else
        somestats.prank = (1+xfitshere)/(nsurrs+1); % rank-based p-value
    end
    if strcmp(leftrightboth,'both')
        somestats.prank = somestats.prank*2; % I think this factor should be in here
    end

    % DO PLOTTING:
%     figure('color','w')
%     plot(statsurr,ones(nsurrs,1),'.k');
%     hold on;
%     plot(statx*ones(2,1),[0,2],'r')
end


end