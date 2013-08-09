% SD_SurrogateTest
% 
% Analyzes the test statistics obtained from surrogate time series compared to
% those measured from the given time series.
% 
% This code was based on information found in:
% "Surrogate data test for nonlinearity including nonmonotonic transforms"
% D. Kugiumtzis Phys. Rev. E 62(1) R25 (2000)
% 
% The generation of surrogates is done by the periphery function,
% SD_MakeSurrogates
% 
% 
% INPUTS:
% x, the input time series
% 
% surrmeth, the method for generating surrogate time series:
%       (i) 'RP': random phase surrogates that maintain linear correlations in
%                 the data but destroy any nonlinear structure through phase
%                 randomization
%       (ii) 'AAFT': the amplitude-adjusted Fourier transform method maintains
%                    linear correlations but destroys nonlinear structure
%                    through phase randomization, yet preserves the approximate
%                    amplitude distribution,
%       (iii) 'TFT': preserves low-frequency phases but randomizes high-frequency phases (as a way of dealing
%                    with non-stationarity, cf.:
%               "A new surrogate data method for nonstationary time series",
%                   D. L. Guarin Lopez et al., arXiv 1008.1804 (2010)
% 
% nsurrs, the number of surrogates to compute (default is 99 for a 0.01
%         significance level 1-sided test)
% 
% extrap, extra parameter, the cut-off frequency for 'TFT'
% 
% teststat, the test statistic to evalute on all surrogates and the original
%           time series. Can specify multiple options in a cell and will return
%           output for each specified test statistic:
%           (i) 'ami': the automutual information at lag 1, cf.
%                 "Testing for nonlinearity in irregular fluctuations with
%                 long-term trends" T. Nakamura and M. Small and Y. Hirata,
%                 Phys. Rev. E 74(2) 026205 (2006)
%           (ii) 'fmmi': the first minimum of the automutual information
%                       function
%           (iii) 'o3': a third-order statistic used in: "Surrogate time
%                 series", T. Schreiber and A. Schmitz, Physica D 142(3-4) 346
%                 (2000)
%           (iv) 'tc3': a time-reversal asymmetry measure. Outputs of the
%                 function include a z-test between the two distributions, and
%                 some comparative rank-based statistics.
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
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = SD_SurrogateTest(x,surrmeth,nsurrs,extrap,teststat)
% Ben Fulcher, 28/1/2011

doplot = 0; % plot outputs to a figure

%% CHECK INPUTS
if nargin < 2 || isempty(surrmeth)
    surrmeth = 'RP'; % randomize phases
end
if nargin < 3 || isempty(nsurrs)
    nsurrs = 99; % create 99 surrogates for a 0.01 significance level 1-sided test
end
if nargin < 4
    extrap = [];
end
if nargin < 5 || isempty(teststat)
    teststat = 'AMI'; % automutual information
end

if ischar(teststat)
    teststat = {teststat};
end

N = length(x); % time-series length

%% OK NOW DO SHIT

% 1) make surrogates
z = SD_MakeSurrogates(x,surrmeth,nsurrs,extrap);
% z is matrix where each column is a surrogate time series

% 2) evaluate test statistic on each surrogate
if ismember('ami1',teststat)
    % look at AMI(1) of surrogates compared to that of signal itself
    % This statistic is used by Nakamura et al. (2006), PRE
    % could use CO_HistogramAMI or TSTL, but I'll use BF_MutualInformation
    % Apparently there are upper and lower bounds on the number of bins to
    % use: [1+log_2(N)], [sqrt(N)]
    nbins = ceil(1+log2(N)); %round(mean([1+log2(N),sqrt(N)]));
    AMIx = BF_MutualInformation(x(1:end-1),x(2:end),'quantile','quantile',nbins);
    AMIsurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        AMIsurr(i) = BF_MutualInformation(z(1:end-1,i),z(2:end,i),'quantile','quantile',nbins);
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
    fmmix = CO_FirstMin(x,'mi');
    fmmisurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        try
            fmmisurr(i) = CO_FirstMin(z(i,:),'mi');
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
    tmp = NL_MS_nlpe(x,de,tau);
    nlpex = tmp.msqerr;
    nlpesurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        res = MS_nlpe(z(:,i),de,tau);
        msqerr = sum(res.^2);
        nlpesurr(i) = msqerr;
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
    tmp = NL_MS_fnn(x,2,1,5,1);
    fnnx = tmp.pfnn_2;
    fnnsurr = zeros(nsurrs,1);
    for i = 1:nsurrs
        tmp = NL_MS_fnn(z(:,i),2,1,5,1);
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
    [~, p, ~, zscore] = ztest(statx, mean(statsurr), std(statsurr), 0.05, leftrightboth);
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
        if (zscstatx < min(xi)) || (zscstatx > max(xi))
            somestats.f = 0; % out of range -- assume p=0 here
        else
            [~, minhere] = min(abs(zscstatx-xi));
            somestats.f = f(minhere); % return probability density where the point is
        end
    end
    
    % What fraction of the range is the sample in?
    medsurr = median(statsurr);
    iqrsurr = iqr(statsurr);
    if iqrsurr == 0
        somestats.mediqr = NaN;
    else
        somestats.mediqr = abs(statx-medsurr)/iqrsurr;
    end

    % rank statistic
    [~, ix] = sort([statx; statsurr]);
    xfitshere = find(ix == 1) - 1;
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
    if doplot
        figure('color','w')
        plot(statsurr,ones(nsurrs,1),'.k');
        hold on;
        plot(statx*ones(2,1),[0,2],'r')
    end
end


end