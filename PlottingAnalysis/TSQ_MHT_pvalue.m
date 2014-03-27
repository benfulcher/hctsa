% ------------------------------------------------------------------------------
% TSQ_MHT_pvalue
% ------------------------------------------------------------------------------
% 
% Returns p-values given a distribution of test stastics and an appropriate null
% distribution. cf. algorithm 18.3, etc. in Hastie's statistical learning text book
% 
%---INPUTS:
%--teststat: the distribution of test statistics from a real labeling of the data
%--teststatR: the distribution of test statistics from shuffling the labels
%             assigned to the data
% 
%---HISTORY:
%--Based on previous code named TSQ_MHT.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function NSig = TSQ_MHT_pvalue(teststat,teststatR,TheTestStatistic,doplot)


% ------------------------------------------------------------------------------
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 3
    TheTestStatistic = 'linclass';
    fprintf(1,'Linear misclassification statistics by default\n');
end
if nargin < 4
    doplot = 1; % plot figures showing results
end

% Test at what type-I error, alpha
alpha = 0.05;

% Set the false discovery rate
fdr = 0.01;

dopool = 1;
if dopool
    fprintf(1,'Pooling test statistics!!!\n');
else
    fprintf(1,'Not pooling test statistics!!!\n')
end

% Variables for use
M = length(teststat); % number of operations independently tested over
K = length(teststatR)/length(teststat); % number of repeats/shuffles in the randomized verison
fprintf(1,['Looks like you used %4.2f shuffles of the data in your randomized version... ' ...
                    'Nice work...\n'],K);

% ------------------------------------------------------------------------------
%% Calculate p-values for shuffled permutations
% ------------------------------------------------------------------------------
pvals = zeros(M,1);
fprintf(1,'Calculating p-values...')

timer = tic;
switch TheTestStatistic
    case 'linclass'
        if dopool
            for i = 1:M
                % pvals(i) = sum(teststat(i) > teststatR(:))/(length(teststatR(:))); % random classifies better
                pvals(i) = sum(teststat(i) > teststatR(:))/(M*K); % random classifies better
                % if (mod(i,floor(M/3))==0)
                %     fprintf(1,'Less than %s remaining! We''re at %u / %u\n',BF_thetime(toc(timer)/M*(M-i)),i,M)
                % end
            end
        else
            for i = 1:M
                pvals(i) = sum(teststat(i)>teststatR(:,i))/K; % random classifies better
            end
        end

    case {'tstat','corrcoef'}
        if dopool
            teststatRcln = teststatR(:);
            for i = 1:M
                pvals(i) = sum(abs(teststat(i))<abs(teststatRcln))/(M*K); % random has more extreme t-statistic
            end
        else
            for i = 1:M
                pvals(i) = sum(abs(teststat(i)) < abs(teststatR(:,i)))/K; % random has more extreme t-statistic
            end
        end
end
fprintf(1,' Done in %s.\n',BF_thetime(toc(timer)))

%% How many are significant?
% Fill structure NSig

% (*) No correction!:
NSig.dumb = sum(pvals < alpha);
fprintf(1,['%u / %u operations significant if we didn''t correct for multiple hypothesis' ...
                ' testing at significance level alpha = %4.2f\n'],NSig.dumb,M,alpha)

% (*) Bonferroni:
NSig.BF = sum(pvals < alpha/M);
fprintf(1,['%u / %u operations significant with Bonferroni correction, M = %u, ' ...
                    'at significance level alpha = %4.2f\n'],NSig.BF,M,M,alpha)

% False discovery rate
[pvals_sort, ix] = sort(pvals,'ascend');
NSig.BH  = find((1:M)'*fdr/M-pvals_sort > 0,1,'last');
if isempty(NSig.BH)
    NSig.BH = 0;
end
fprintf(1,['%u / %u operations significant' ...
        ' with False Discovery Rate, FDR at most = %4.2f' ...
        ' using the Benjamini-Hochberg method\n'],NSig.BH,M,fdr)
if NSig.BH > 0
    fprintf(1,['Correponds to a threshold of p = %4.2f and a test statistic' ...
                    ' %4.2f\n'],pvals_sort(NSig.BH),teststat(ix(NSig.BH)))
end

if doplot
    % False discovery rate plot
    figure('color','w'); box('on'); hold on
    cc = BF_getcmap('set1',3,1);
    plot(fdr.*(1:M)/M,'-','color',cc{3},'linewidth',2);
    plot(pvals_sort(1:NSig.BH),'o','MarkerFaceColor',cc{1});
    plot((NSig.BH+1:length(pvals)),pvals_sort(NSig.BH+1:end),'o','MarkerFaceColor',cc{2});
    xlabel('Operations ordered by p-value')
    ylabel('p-value');
    ylim([0,fdr])
    legend({'linear fdr threshold',['Significant: FDR<' num2str(fdr)],['Not significant: FDR > ' num2str(fdr)]})
end

%% Estimate False discovery rate
% (algorithm 18.3 in Hastie)
% choose cut-point, C
if (NSig.BH > 0) % some were significant
    C = teststat(ix(NSig.BH)); % classification error less than the BH threshold
else
    fprintf(1,'None significant... ok...\n');
    C = teststat(ix(1));
end

switch TheTestStatistic
    case 'linclass'
        Robs = sum(teststat <= C); % we observe significant
        EV = sum(teststatR(:) <= C)/K; % randomly significant
    case {'tstat','corrcoef'}
        Robs = sum(abs(teststat) >= abs(C)); % we observe significant
        EV = sum(abs(teststatR(:)) >= abs(C))/K; % randomly significant
end
FDR = EV/Robs;
fprintf(1,'False discovery rate at threshold C = %4.2f is %4.2f\n',C,FDR)

%% Plot histogram with significance thresholds included
if doplot
    plothow = 'threshold';
    nbins = 25; % number of bins for histogram

    figure('color','w'); box('on'); hold on
    switch TheTestStatistic
        case 'linclass'
            allstats = [teststat;teststatR(:)];
            x = linspace(min(allstats),max(allstats),nbins+1);
        case 'tstat'
            x = linspace(-1,1,nbins+1);
%             Maxmax = max([max(abs(min(teststatR(:))),max(teststatR(:))),max(abs(min(teststat)),max(teststat))]);
%             x = linspace(-Maxmax,Maxmax,nbins+1);
        case 'corrcoef'
            mmm = max(abs(teststat));
            if mmm > 0.7
                mmm = 1; % span histograms from -1,1
            end
            x = linspace(-mmm,mmm,nbins+1);
    end
    
    x = x(1:end-1)+(x(2)-x(1))/2; % centre bins
    NR = hist(teststatR(~isnan(teststatR)),x)/K; % randomized
    switch plothow
        case 'threshold'
            lw = 2;
            colormap(BF_getcmap('set1',2,0))
            N = hist(teststat(~isnan(teststat)),x);
            MaxMax = max([N,NR]);
            bar(x,[N;NR]','grouped');
            cc = BF_getcmap('set2',4,1);
            
            % False Discovery Rates
            fdrr = [0.05,0.01,0.001];
            nsigBHs = zeros(length(fdrr),1);
            for i = 1:length(fdrr)
                YepHere = find((1:M)'*fdrr(i) / M-pvals_sort > 0,1,'last');
                if isempty(YepHere)
                    nsigBHs(i) = 0;
                else
                    nsigBHs(i) = YepHere;
                end
                if (nsigBHs(i) > 0) % some significant
                    Ci = teststat(ix(nsigBHs(i)));
                    switch TheTestStatistic
                        case {'tstat','corrcoef'}
                            plot([Ci,Ci],[0,MaxMax],'--','color',cc{1+i},'LineWidth',lw)
                            plot(-[Ci,Ci],[0,MaxMax],'--','color',cc{1+i},'LineWidth',lw)
    %                         text(Ci,MaxMax*(1-0.1*(i+1)),['FDR=' num2str(fdrr(i))]);
                        case 'linclass'
                            plot([Ci,Ci],[0,MaxMax],'--','color',cc{1+i},'LineWidth',lw)
    %                         text(Ci,MaxMax*(1-0.1*(i+1)),['FDR=' num2str(fdrr(i))]);
                    end
                end
                
            end
            
            % Bonferroni
            CBF = teststat(ix(find(pvals_sort < alpha/M,1,'last')));
            if ~isempty(CBF)
                switch TheTestStatistic
                    case {'tstat','corrcoef'}
                        plot([CBF,CBF],[0,MaxMax],'--','color',cc{1},'LineWidth',lw)
                        plot(-[CBF,CBF],[0,MaxMax],'--','color',cc{1},'LineWidth',lw)
                    case 'linclass'
                        plot([CBF,CBF],[0,MaxMax],'--','color',cc{1},'LineWidth',lw)
            %             text(CBF,MaxMax,'Bonferroni');
                end
            end
            
            % Legend
            switch TheTestStatistic
                case {'tstat','corrcoef'}
                    legend(['Unshuffled (' num2str(M) ')'],...
                    [num2str(K) ' x shuffled (pooled)'],...
                    ['FDR <= ' num2str(fdrr(1)) ' (' num2str(nsigBHs(1)) ')'],'', ...
                    ['FDR <= ' num2str(fdrr(2)) ' (' num2str(nsigBHs(2)) ')'],'', ...
                    ['FDR <= ' num2str(fdrr(3)) ' (' num2str(nsigBHs(3)) ')'],'', ...
                    ['Bonferroni \alpha = ' num2str(alpha) ' (' num2str(NSig.BF) ')'],'');
                case 'linclass'
                    legend(['Unshuffled (' num2str(M) ')'],...
                    [num2str(K) ' x shuffled (pooled)'],...
                    ['FDR <= ' num2str(fdrr(1)) ' (' num2str(nsigBHs(1)) ')'], ...
                    ['FDR <= ' num2str(fdrr(2)) ' (' num2str(nsigBHs(2)) ')'], ...
                    ['FDR <= ' num2str(fdrr(3)) ' (' num2str(nsigBHs(3)) ')'], ...
                    ['Bonferroni \alpha = ' num2str(alpha) ' (' num2str(NSig.BF) ')']);
            end
            
        case '3classes'
            colormap(BF_getcmap('set1',3,0))
            % plot significant as its own histogram
            switch TheTestStatistic
                case {'tstat','corrcoef'}
                    Ns = hist(teststat(~isnan(teststat(abs(teststat)<abs(C)))),x);
                    Nns = hist(teststat(~isnan(teststat(abs(teststat)>=abs(C)))),x);
                case 'linclass'
                    Ns = hist(teststat(~isnan(teststat(teststat<C))),x);
                    Nns = hist(teststat(~isnan(teststat(teststat>C))),x);
            end
            bar(x,[Nns;NR;Ns]','grouped');%,'FaceColor','w','EdgeColor','k');
            legend(['Not significant (' num2str(M-NSig.BH) ')'],...
                [num2str(K) ' x shuffled (pooled)'],...
                ['Significant, FDR <= ' num2str(fdr) ' (' num2str(NSig.BH) ')'])
                
        case '3colors'
            colormap(BF_getcmap('set1',3,0))
            N1 = hist(teststat(~isnan(teststat)),x); %/sum(~isnan(teststat));
            % plot significant as a different class
            N1ns = N1; N1s = N1;
            switch TheTestStatistic
                case {'tstat','corrcoef'}
                    N1ns(abs(x) >= abs(C)) = 0; % not significant
                    N1s(abs(x) < abs(C)) = 0; % significant
                case 'linclass'
                    N1ns(x<C) = 0; % not significant
                    N1s(x>C) = 0; % significant
            end
            bar(x,[N1ns; NR; N1s]','grouped');%,'FaceColor','w','EdgeColor','k');
            legend(['Not significant (' num2str(M-NSig.BH) ')'],...
                [num2str(K) ' x shuffled (pooled)'],...
                ['Significant, FDR <= ' num2str(fdr) ' (' num2str(NSig.BH) ')'])
    end
    
    set(gca,'ylim',[0,MaxMax])
    ylabel('Frequency')
    set(gcf,'Position',[400,600,500,270]);
    switch TheTestStatistic
        case 'linclass'
            xlabel('Linear Classification Rate (%)')
        case 'tstat'
            xlabel('t statistic')
        case 'corrcoef'
            xlabel('Linear Correlation Coefficient, R')
            set(gca,'xlim',[-mmm,mmm])
    end
    mywd = pwd;
    islash = regexp(mywd,'/');
    title(mywd(islash(end)+1:end),'interpreter','none');
end

end