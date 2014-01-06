% TSQ_IndividualFeatures
% 
% Searches for individual features that help distinguish a known classification
% of the time series
%
%--INPUTS
%-ClassMeth: the classification method to use. The following are possible options:
%            'linclasscv', 'linclass', 'linclasscv', 'knn', 'knncv'
%
%
%--HISTORY
%% uses BioInformatics toolbox functions
% Ben Fulcher 13/9/2010
% Ben Fulcher 28/9/2010 added subset input
% Ben Fulcher 18/3/2011 changed reverse to cvalid -- specify cross-validation
% options as in TSQ_cfnerr: e.g, {'kfold',10,2}; 'ttest' (default) -- Absolute
% value two-sample t-test with pooled variance estimate. 'entropy' -- Relative
% entropy, also known as Kullback-Leibler distance or divergence.
% 'bhattacharyya' -- Minimum attainable classification error or Chernoff bound.
% 'roc' -- Area between the empirical receiver operating characteristic (ROC)
% curve and the random classifier slope.
% 'wilcoxon' -- Absolute value of the
% u-statistic of a two-sample unpaired Wilcoxon test, also known as
% Mann-Whitney.
% Rperm is the permutation used if randomize is specified: no
% longer an output
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

function [ifeat, teststat, testspread] = TSQ_IndividualFeatures(WhatData,ClassMethod,cvalid,randomize,plotoutputs)

%% Check inputs

if nargin < 1 || isempty(WhatData)
    WhatData = 'norm';
end
if nargin < 2 || isempty(ClassMethod)
    ClassMethod = 'linclass';
    fprintf(1,'Using ''%s'' by default\n', ClassMethod);
end
if nargin < 3
    cvalid = [];
end
if nargin < 4
    randomize = 0;
end
if nargin < 5
    plotoutputs = 5; % number of figures outputs by default
end

% set plotting options in a structure: plotopts
if ~isstruct(plotoutputs)
    if plotoutputs > 0 % plot all outputs
        plotopts.histogram = 1;
        plotopts.nfigs = plotoutputs;
    else
        plotopts.histogram = 0;
        plotopts.nfigs = 0;
    end
else
    if isfield(plotoutputs,'histogram')
        plotopts.histogram = plotoutputs.histogram; % plot histogram
    else
        plotopts.histogram = 1;
    end
    if isfield(plotoutputs,'nfigs')
        plotopts.nfigs = plotoutputs.nfigs;
    else
        plotopts.nfigs = 5;
    end
end

%% Load the data
% if specified a string: 'norm' or 'cl' will retrieve
% otherwise inputted the data matrix as TS_DataMat and given manual groupings in 
% kwgs and in gi
switch WhatData
case 'norm'
    TheDataFile = 'HCTSA_N.mat';
case 'cl'
    TheDataFile = 'HCTSA_cl.mat';
end
fprintf(1,'Loading data from %s...',TheDataFile);
load(TheDataFile,'TS_DataMat','Operations','TimeSeries');
fprintf(1,' Loaded.\n');

% Randomize the data matrix to check null
if randomize % shuffle elements of the data matrix
    fprintf(1,'Randomly permuting the group information as a null comparison...\n');
    Rperm = randperm(size(TS_DataMat,1)); % this is an output
    TS_DataMat = TS_DataMat(Rperm,:);
end

%% Run the algorithm
gig = [TimeSeries.Group]; % Use group form
switch ClassMethod
    case 'linclasscv' % quantify linear classification for each
        teststat = zeros(size(TS_DataMat,2),1); % cross-validation classification rates
        testspread = zeros(size(TS_DataMat,2),1); % standard deviation in cross-validation classification rates
        
        if isempty(cvalid) % set default cross-validation
            cvalid = {'kfold',10,1}; % 10-fold CV with no repeats
        end
%         cp = cvpartition(gig,'kfold',cvalid{2}); % ?-fold stratified cross-validation
%         cp = cvpartition(gig,'k',10); % 10-fold stratified cross-validation

        nrepeats = cvalid{3}; % make this many different 10-fold partitions of the data
                              % ensures they're the same for all operations
        cps = cell(nrepeats,1);
        for i = 1:nrepeats;
            cps{i} = cvpartition(gig,'kfold',cvalid{2}); % ?-fold stratified cross-validation;
        end
        
        fprintf(1,['DOING CROSS VALIDATION WITH %u REPEATS USING ' ...
                        '%u-fold CV'],nrepeats,cvalid{2});
        
        fn_classify = @(XT,yT,Xt,yt) sum(yt~=classify(Xt,XT,yT,'linear'))/length(yt);
        times = zeros(size(TS_DataMat,2),1); % time each run:
        for i = 1:size(TS_DataMat,2)
            tic
            try
%                 errs = TSQ_cfnerr('classify','linear',TS_DataMat(:,i),gig,[],cvalid);
                te = cell(nrepeats,1);
                for j = 1:nrepeats
                    % use partition cps{j}
                    te{j} = crossval(fn_classify,TS_DataMat(:,i),gig,'partition',cps{j});
%                     errs = TSQ_cfnerr('classify','linear',TS_DataMat(:,i),gig,[],cp,cvalid{3});
                end
                errs = vertcat(te{:}); % agglomerate across repeats
                teststat(i) = mean(errs);
                testspread(i) = std(errs);
%                 teststat(i) = crossval('mcr',TS_DataMat(:,i),gig,'predfun',classf,'partition',cp);
%                 mcrs = crossval(F_linclass,TS_DataMat(:,i),gig,'partition',cp);
%                 This code with 'mcr' is the same as mean(mcrs)
            catch emsg
                disp(emsg)
                teststat(i) = NaN;
            end
            times(i) = toc;
            if mod(i,floor(size(TS_DataMat,2)/4))==0
                disp(['Less than ' benrighttime(mean(times(1:i))*(size(TS_DataMat,2)-i)) ...
                    ' remaining! We''re at ' num2str(i) ' / ' num2str(size(TS_DataMat,2))])
            end
        end
        [teststat, ifeat] = sort(teststat,'ascend');
        testspread = testspread(ifeat)*100; % sort, convert to percentages
        teststat = teststat*100; % convert to percentages
        
    case 'linclass'
        fprintf(1,['Comparing the performance of %u operations using ' ...
                'in-sample linear classification WITHOUT cross validation...\n'],length(Operations))
        timer = tic;
        teststat = zeros(size(TS_DataMat,2),1); % in-sample classification rates
        for i = 1:size(TS_DataMat,2)
            try
                [~,err] = classify(TS_DataMat(:,i),TS_DataMat(:,i),gig,'linear'); % in-sample errors
                teststat(i) = err;
            catch
                teststat(i) = NaN;
            end
        end
        fprintf(1,'Done. Took %s.\n',BF_thetime(toc(timer)));
        [teststat,ifeat] = sort(teststat,'ascend');
        teststat = teststat*100;
        testspread = [];
        
    case {'knn','knn_matlab'}
        k = 3;
        fprintf(1,'IN-SAMPLE BEN KNN(%u)\n',k)
        teststat = zeros(size(TS_DataMat,2),1); % in-sample classification rates
        for i = 1:size(TS_DataMat,2) % for each feature individually
            [~,err] = benknn(TS_DataMat(:,i),TS_DataMat(:,i),3,gig,gig); % classifies in-sample
            teststat(i) = err;
        end
        [teststat,ifeat] = sort(teststat,'ascend');
        teststat = teststat*100;
        testspread = [];
        
    case 'knncv'
        teststat = zeros(size(TS_DataMat,2),1); % classification rates
        testspread = zeros(size(TS_DataMat,2),1); % standard deviation in cross-validation classification rates
        
        k = 3; % number of nearest neighbors
        if isempty(cvalid) % set default cross-validation
            cvalid = {'kfold',10,1}; % 10-fold CV with no repeats
        end
        
%         cp = cvpartition(gig,'kfold',nfolds); % statified k-fold crossvalidation
        disp([num2str(cvalid{1}) ' ' cvalid{1} ' cross validation BEN KNN(' num2str(k) ')'])
        for i = 1:size(TS_DataMat,2) % for each feature individually
            err = TSQ_cfnerr('knn',k,TS_DataMat(:,i),gig,[],cvalid);
            teststat(i) = mean(err);
            testspread(i) = std(err);
%             errs = zeros(nfolds,1);
%             for j = 1:nfolds
%                 itrain = find(training(cp,j)); % indicies for training
%                 itest = find(test(cp,j)); % indicies for testing
%                 [~,errs(j)] = benknn(TS_DataMat(itrain,i),TS_DataMat(itest,i),k,gig(itrain),gig(itest)); % classifies test data
%             end
%             teststat(i) = mean(errs);
        end
        [teststat,ifeat] = sort(teststat,'ascend');
        testspread = testspread(ifeat)*100; % sort, convert to percentages
        teststat = teststat*100; % convert to percentages
        
    case {'knn_spider','svm'}
        % define the spider model
        if strcmp(ClassMethod,'knn')
            a = SPIDER_getmemodel('knn',3); % define a knn(3) model
            disp('using knn(3)')
        else
            a = SPIDER_getmemodel('svm',{{'linear'}}); % define a svm (linear) model
            disp('using svm(linear)')
        end
        
        % initialize the variables
%         ifeat = zeros(nbest,1); % stores indicies of features chosen at each stage
%         teststat = cell(nbest,1); % cross-validation classification rates for all features at each stage
        
        % calculate losses across all features (in combination with those
        %                                   already chosen)
        teststat = zeros(size(TS_DataMat,2),1); % in-sample classification rates

        for i = 1:size(TS_DataMat,2)        
            dtrain = makeitdata(TS_DataMat(:,i),gi);
            try
                [~,a] = train(a,dtrain); % train classification model on training data
                lossme = loss(test(a,dtrain),'class_loss'); % get in-sample loss
                teststat(i) = lossme.Y;
            catch
                teststat(i) = NaN;
            end
        end
        [teststat,ifeat] = sort(teststat,'ascend');
        teststat = teststat*100;
        
    otherwise
        fprintf(1,'Hello!\n')
        [ifeat,teststat] = rankfeatures(TS_DataMat',gig,'criterion',ClassMethod);
        [teststat,ix] = sort(teststat,'descend');
        ifeat = ifeat(ix);
        testspread = [];
end
    
% Print the top 25
topn = min(25,length(Operations));
for i = 1:topn
    fprintf(1,['[%u] {%u} %s -- %s :: %4.2f%%\n'],ifeat(i),Operations(ifeat(i)).ID, ...
                    Operations(ifeat(i)).Name,Operations(ifeat(i)).Keywords,teststat(i));
end



%-----------------------This function no longer plots outputs-----------------
% BF, 2014-01-06

%% Plot outputs
if plotopts.histogram
    % 1) a figure to show the distribution of test statistics across all
    % features:
    % ksdensity
%     figure('color','w'); box('on'); hold on
%     [f,xi] = ksdensity(teststat);
%     plot(xi,f,'b');
%     plot(teststat(1),0,'or');
%     xlabel('individual classification errors')
%     ylabel('probability density')
    % histogram
    figure('color','w');
    hist(teststat,10);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    % bar(xout,n)
    xlabel('individual classifiction error')
    ylabel('frequency')
end
% 
% if plotopts.nfigs > 0
%     ntoplot = 5; % subplots per figure
%     howmanyrepeats = plotopts.nfigs; % number of plots to produce
%     
%     % classfun = @(XT,yT,Xt,yt) (sum(~strcmp(yt,classify(Xt,XT,yT,'linear'))));
%     % cvMCR = crossval('mcr',X,y,'predfun',fun);
%     labels.op_ids = m_ids;
%     labels.mlab = mlab;
%     labels.mkw = mkw;
%     for k = 1:howmanyrepeats
%         if k*ntoplot>length(ifeat)
%             break % we've exceeded number of features
%         end
%         r = ((k-1)*ntoplot+1:k*ntoplot);
%         feathere = ifeat(r);
%         labels.teststat = teststat(r);
%         TSQ_plot_kss(kwgs,gi,TS_DataMat,feathere,labels,k==1);
%         % plot kernel smoothed distributions
%     end
% 
% %     figure('color','w')
% %     for i = 1:ntoplot
% %         subplot(ntoplot,1,i); hold on; box('on')
% %         i = ntoplot+i; % plot also the next best
% %         for j=1:Ng
% %             plot_ks(TS_DataMat(gi{j},ifeat(i)),c{j},0) % plot the jth group
% %         end
% % 
% %     %     % get cross-validation classification rate
% %     %     cp = cvpartition(gig,'k',10); % 10-fold stratified cross-validation
% %     %     cvMCR = crossval('mcr',TS_DataMat(:,ifeat(i)),gig,'predfun',classf,'partition',cp);
% % 
% %         title(['[' num2str(m_ids(ifeat(i))) ']' mlab{ifeat(i)} ' [' mkw{ifeat(i)} '] -- ' ClassMethod ' = ' num2str(teststat(i))],'interpreter','none')
% %         set(gca,'YTick',[])
% %     %     if i < i+ ntoplot
% %     %         set(gca,'XTick',[]);
% %     %     end
% %         if i == ntoplot + 1
% %             leglabs = cell(Ng*2,1);
% %             for j = 1:Ng
% %                 leglabs{2*(j-1)+1} = [kwgs{j} ' (' num2str(length(gi{j})) ')'];
% %                 leglabs{2*(j-1)+2} = '';
% %             end
% %             legend(leglabs)
% %         end
% %     end
% end

% function plot_ks(v,c,swap)
%     % vector v is the vector of a given group
%     % c is the color
%     if nargin<3, swap = 0; end
%     [f x] = ksdensity(v,linspace(min(v),max(v),1000),'function','pdf');
% %     [f x] = ksdensity(v,'function','pdf');
%     r = zeros(length(v),1);
%     for m=1:length(v);
%         r(m) = find(x>=v(m),1,'first');
%     end
%     % crop to just 50 points across the range:
% %     r = sort(r,'ascend'); r=r(round(linspace(1,length(r),50)));
%     if swap
%         plot(-f,x,'color',c);
%         plot(-f(r),x(r),'.','color',c)
%     else
%         plot(x,f,'color',c);
%         plot(x(r),f(r),'.','color',c)
%     end
% end


end