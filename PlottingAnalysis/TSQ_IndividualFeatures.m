% ------------------------------------------------------------------------------
% TSQ_IndividualFeatures
% ------------------------------------------------------------------------------
% 
% Searches for individual features that help distinguish a known classification
% of the time series
%
%---INPUTS
% ClassMeth: the classification method to use. The following are possible options:
%            'linclasscv', 'linclass', 'linclasscv', 'knn', 'knncv'
%
%---HISTORY
% Previously TSQ_ML_rankfeatures
% uses BioInformatics toolbox functions
% Ben Fulcher 13/9/2010
% Ben Fulcher 28/9/2010 added subset input
% Ben Fulcher 18/3/2011 changed reverse to CrossVal -- specify cross-validation
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

function [ifeat, teststat, testspread] = TSQ_IndividualFeatures(WhatData,ClassMethod,CrossVal,randomize,plotoutputs)

% --------------------------------------------------------------------------
%%                          Check inputs
% --------------------------------------------------------------------------

if nargin < 1 || isempty(WhatData)
    WhatData = 'norm';
end
if nargin < 2 || isempty(ClassMethod)
    ClassMethod = 'linclass';
    fprintf(1,'Using ''%s'' by default\n', ClassMethod);
end
if nargin < 3
    CrossVal = {'kfold',10,1};
    % CrossVal = {'leaveout'};
end
if nargin < 4
    randomize = 0;
end
if nargin < 5
    plotoutputs = 5; % number of figures outputs by default
end

% --------------------------------------------------------------------------
% Place plotting options in a structure: plotopts
% --------------------------------------------------------------------------
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

% --------------------------------------------------------------------------
%%                          Load the data
% --------------------------------------------------------------------------
if ischar(WhatData)
    switch WhatData
    case 'norm'
        TheDataFile = 'HCTSA_N.mat';
    case 'cl'
        TheDataFile = 'HCTSA_cl.mat';
    end
    fprintf(1,'Loading data from %s...',TheDataFile);
    load(TheDataFile,'TS_DataMat','Operations','TimeSeries');
    fprintf(1,' Loaded.\n');
elseif isstruct(WhatData)
    % Already loaded and given here as a structure
    TS_DataMat = WhatData.TS_DataMat;
    Operations = WhatData.Operations;
    TimeSeries = WhatData.TimeSeries;
    fprintf(1,'Data provided, adapted successfully.\n');
end


% --------------------------------------------------------------------------
%      Randomize the rows of the data matrix for a null distribution
% --------------------------------------------------------------------------
if randomize % shuffle elements of the data matrix
    fprintf(1,'Randomly permuting the group information as a null comparison...\n');
    Rperm = randperm(length(TimeSeries)); % this is an output
    TsGroups = [TimeSeries.Group];
    % This ([TimeSeries.Group] = TsGroups(Rperm);) doesn't work, so I have to be really inefficient:
    for i = 1:length(TimeSeries)
        TimeSeries(i).Group = TsGroups(Rperm(i)); % Rewrite with random permutation
    end
    % Another (equivalent) option to avoid this loop is to permute the data matrix...
    % Rperm = randperm(size(TS_DataMat,1)); % this is an output
    % TS_DataMat = TS_DataMat(Rperm,:);
end


% --------------------------------------------------------------------------
%% Define the train/test classification rate functions
% --------------------------------------------------------------------------
switch ClassMethod
case {'linear','linclass'}
    fprintf(1,'A linear classifier\n');
    fn_classify = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest == classify(Xtest,XTrain,yTrain,'linear'))/length(ytest);
case 'diaglinear'
    fprintf(1,'A Naive Bayes classifier\n');
    fn_classify = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest == classify(Xtest,XTrain,yTrain,'diaglinear'))/length(ytest);
case {'svm','svmlinear'}
    fprintf(1,'A linear support vector machine\n');
    fn_classify = @(XTrain,yTrain,Xtest,ytest) ...
                    sum(ytest == svmclassify(svmtrain(XTrain,yTrain, ...
                                'Kernel_Function','linear'),Xtest))/length(ytest);
otherwise
    error('Unknown classification method ''%s''',ClassMethod)
end

% --------------------------------------------------------------------------
%% Cross-validation data partition functions
% --------------------------------------------------------------------------

TimeSeriesGroup = [TimeSeries.Group]'; % Use group form

switch CrossVal{1}
case 'kfold'
    Nfolds = CrossVal{2};
    nrepeats = CrossVal{3}; % make this many different 10-fold partitions of the data
                          % ensures they're the same for all operations

    fprintf(1,['Doing cross validation with %u repeats using ' ...
                  '%u-fold cross validation...\n'],nrepeats,Nfolds);

    % ?-fold stratified cross-validation;
    fn_partition = @(TheLabels) cvpartition(TheLabels,'kfold',Nfolds);

    % Determine partitions ahead of time first, so you do the same set of partitions
    % at each iteration (for consistency):
    Partitions = cell(nrepeats,1);
    for i = 1:nrepeats;
        Partitions{i} = fn_partition(TimeSeriesGroup); 
    end
    
case 'leaveout'
    fn_partition = @(TheLabels) cvpartition(TheLabels,'leaveout');
    fprintf(1,'Doing leave-one-out cross validation\n');
end


% --------------------------------------------------------------------------
%%                     Loop over all features
% --------------------------------------------------------------------------

switch CrossVal{1}
case 'kfold'
    MeanClassificationRate = zeros(size(TS_DataMat,2),1);
    timer = tic;
    for i = 1:size(TS_DataMat,2)
        TestRate = zeros(nrepeats,1);
        for j = 1:nrepeats
            TestRate(j) = mean(crossval(fn_classify,TS_DataMat(:,i), ...
                                        TimeSeriesGroup, ...
                                        'partition',Partitions{j}));
            % Mean across k folds
        end
        MeanClassificationRate(i) = mean(TestRate); % Mean of means (across all repartitions) 
        
        if (mod(i,floor(size(TS_DataMat,2)/4))==0)
            fprintf(1,'Less than %s remaining! We''re at %u / %u\n', ...
                        BF_thetime(toc(timer)/i*(size(TS_DataMat,2)-i)),i,size(TS_DataMat,2))
        end
    end

    teststat = MeanClassificationRate*100; % Convert to percentages

case 'leaveout'
    MeanClassificationRate = zeros(size(TS_DataMat,2),1);
    timer = tic;
    for i = 1:size(TS_DataMat,2)
        try
            MeanClassificationRate(i) = mean(crossval(fn_classify,TS_DataMat(:,i), ...
                                TimeSeriesGroup,'partition',fn_partition(TimeSeriesGroup)));
        catch
            MeanClassificationRate(i) = NaN;
        end
        
        if (mod(i,floor(size(TS_DataMat,2)/4))==0)
            fprintf(1,'Less than %s remaining! We''re at %u / %u\n', ...
                        BF_thetime(toc(timer)/i*(size(TS_DataMat,2)-i)),i,size(TS_DataMat,2))
        end
    end
    teststat = MeanClassificationRate*100; % Convert to percentages
    
case 'none' % in-sample
    fprintf(1,'Comparing %u operations using in-sample classification...',length(Operations))
    
    timer = tic;
    teststat = zeros(size(TS_DataMat,2),1); % in-sample misclassification rates
    for i = 1:size(TS_DataMat,2)
        try
            % Use the same data for training and testing:
            teststat(i) = fn_classify(TS_DataMat(:,i),TimeSeriesGroup,TS_DataMat(:,i),TimeSeriesGroup);
        catch
            teststat(i) = NaN;
        end
    end
    fprintf(1,' Done in %s.\n',BF_thetime(toc(timer)));
    
    % [teststat, ifeat] = sort(teststat,'descend');
    teststat = teststat*100; % Convert to percentages
    % testspread = []; % no spread because no cross-validation/repeats
    
    % Give mean and that expected from random classifier (maybe a little overfitting)
    fprintf(1,'Mean across %u operations = %4.2f; (Random guessing for %u equiprobable classes = %4.2f)\n', ...
            length(Operations),mean(teststat),length(unique(TimeSeriesGroup)),100/length(unique(TimeSeriesGroup)));

%     case {'knn','knn_matlab'}
%         k = 3;
%         fprintf(1,'IN-SAMPLE BEN KNN(%u)\n',k)
%         teststat = zeros(size(TS_DataMat,2),1); % in-sample classification rates
%         for i = 1:size(TS_DataMat,2) % for each feature individually
%             [~,err] = benknn(TS_DataMat(:,i),TS_DataMat(:,i),3,TimeSeriesGroup,TimeSeriesGroup); % classifies in-sample
%             teststat(i) = err;
%         end
%         [teststat,ifeat] = sort(teststat,'ascend');
%         teststat = teststat*100;
%         testspread = [];
%         
%     case 'knncv'
%         teststat = zeros(size(TS_DataMat,2),1); % classification rates
%         testspread = zeros(size(TS_DataMat,2),1); % standard deviation in cross-validation classification rates
%         
%         k = 3; % number of nearest neighbors
%         if isempty(CrossVal) % set default cross-validation
%             CrossVal = {'kfold',10,1}; % 10-fold CV with no repeats
%         end
%         
% %         cp = cvpartition(TimeSeriesGroup,'kfold',nfolds); % statified k-fold crossvalidation
%         fprintf(1,'%s cross validation BEN KNN(%u)\n',CrossVal{1},k)
%         for i = 1:size(TS_DataMat,2) % for each feature individually
%             err = TSQ_cfnerr('knn',k,TS_DataMat(:,i),TimeSeriesGroup,[],CrossVal);
%             teststat(i) = mean(err);
%             testspread(i) = std(err);
%         end
%         [teststat,ifeat] = sort(teststat,'ascend');
%         testspread = testspread(ifeat)*100; % sort, convert to percentages
%         teststat = teststat*100; % convert to percentages
    %     
    % otherwise
    %     fprintf(1,'Hello!\n')
    %     [ifeat, teststat] = rankfeatures(TS_DataMat',TimeSeriesGroup,'criterion',ClassMethod);
    %     [teststat, ix] = sort(teststat,'descend');
    %     ifeat = ifeat(ix);
    %     testspread = [];
end

% --------------------------------------------------------------------------
%%          Display information the top n operations
% --------------------------------------------------------------------------
[teststat_sort, ifeat] = sort(teststat,'descend');

topn = min(10,length(Operations));
for i = 1:topn
    fprintf(1,'[%u] {%u} %s -- %s :: %4.2f%%\n',ifeat(i),Operations(ifeat(i)).ID, ...
            Operations(ifeat(i)).Name,Operations(ifeat(i)).Keywords,teststat_sort(i));
end


% --------------------------------------------------------------------------
%%                          Plot outputs
% --------------------------------------------------------------------------
if plotopts.histogram
    % 1) a figure to show the distribution of test statistics across all
    % features:

    figure('color','w');
    hist(teststat,10);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    % bar(xout,n)
    xlabel('Individual classification error')
    ylabel('Frequency')
end

if plotopts.nfigs > 0
    ntoplot = 5; % subplots per figure
    howmanyrepeats = plotopts.nfigs; % number of plots to produce

    labels.op_ids = m_ids;
    labels.mlab = mlab;
    labels.mkw = mkw;
    for k = 1:howmanyrepeats
        if k*ntoplot > length(ifeat)
            break % We've exceeded number of features
        end
        r = ((k-1)*ntoplot+1:k*ntoplot);
        feathere = ifeat(r);
        labels.teststat = teststat(r);
        
        % Go through and plot
        figure('color','w');
        for Opi = 1:TopHowMany
            BestHere = ix(Opi);
            if strcmp(My_Region,'all')
                Region_ind = TheRegions(BestHere);
                Op_ind = TheOperations(BestHere);
            else
                Region_ind = My_Region;
                Op_ind = BestHere;
            end

            % Define the values of the full feature vector:
            FullFV = TS_DataMat(([TimeSeries.RegionID] == Region_ind),Op_ind);
    
            fprintf(1,['Best classification rate was %4.2f%% for ' ...
                            'feature %s in region %u\n'],ClassificationRates(Region_ind,Op_ind), ...
                            Operations(Op_ind).Name,Region_ind);
            fprintf(1,'Has a p-value of %5.5f\n',pvalue_pooled(Region_ind,Op_ind));

            subplot(TopHowMany,1,Opi)
            box('on');

            switch ks_or_hist
            case 'ks'
                hold on
                for i = 1:2
                    TheFeatureVectorClassi = TS_DataMat(([TimeSeries.Group]==i) &  ...
                                        ([TimeSeries.RegionID] == Region_ind),Op_ind);
                    [f, x] = ksdensity(TheFeatureVectorClassi);
                    plot(x,f,'color',MarkerColors{i},'LineWidth',2);
                    % Add dots:
                    r = (arrayfun(@(m)find(x >= m,1,'first'),TheFeatureVectorClassi));
                    plot(x(r),f(r),'o','MarkerFaceColor',MarkerColors{i},'MarkerEdgeColor',MarkerColors{i})
                end
                legend('Control','','Schiz','')
            case 'hist'
                nbins = 10;
                TheFeatureVectorClass1 = TS_DataMat(([TimeSeries.Group]==1) &  ...
                                        ([TimeSeries.RegionID] == Region_ind),Op_ind);
                TheFeatureVectorClass2 = TS_DataMat(([TimeSeries.Group]==2) &  ...
                                        ([TimeSeries.RegionID] == Region_ind),Op_ind);

                x = linspace(min(FullFV),max(FullFV),nbins+1);
                x = x(1:end-1)+(x(2)-x(1))/2; % centre bins
                N1 = hist(TheFeatureVectorClass1(~isnan(TheFeatureVectorClass1)),x);
                N2 = hist(TheFeatureVectorClass2(~isnan(TheFeatureVectorClass2)),x);
                histwidth = x(2)-x(1);
                if histwidth > 0
                    plotbar = @(z,N,c) rectangle('Position',[z,0,histwidth/2,N],'FaceColor',c);
                    for k = 1:length(x);
                        if N1(k)>0
                            plotbar(x(k)-histwidth/2,N1(k),MarkerColors{1});
                        end
                        if N2(k)>0
                            plotbar(x(k),N2(k),MarkerColors{2});
                        end
                    end
                    colormap(vertcat(MarkerColors{:}));
                    legend('Control','Schiz')
                end
            end
    
            % Set xlim:
            if min(FullFV) < max(FullFV)
                xlim([min(FullFV),max(FullFV)])
            end


            xlabel(sprintf('Cluster of %u centered at %s', ...
                    length(opid_cl{ismember(cellfun(@(x)x(1),opid_cl),[Operations(Op_ind).ID])}), ...
                                    Operations(Op_ind).Name),'interpreter','none')
            ylabel('Probability density')
            title(sprintf(['%4.2f%% classification (p = %5.5f) in region %u\n'], ...
                                    ClassificationRates(Region_ind,Op_ind), ...
                                    pvalue_pooled(Region_ind,Op_ind),Region_ind), ...
                                    'interpreter','none');
            set(gcf,'Position',[987,98,290,1047])
        end
        
        
        
        
        
        % TSQ_plot_kss(kwgs,gi,TS_DataMat,feathere,labels,k==1);
        % plot kernel smoothed distributions
    end
end
%     figure('color','w')
%     for i = 1:ntoplot
%         subplot(ntoplot,1,i); hold on; box('on')
%         i = ntoplot+i; % plot also the next best
%         for j=1:Ng
%             plot_ks(TS_DataMat(gi{j},ifeat(i)),c{j},0) % plot the jth group
%         end
% 
%     %     % get cross-validation classification rate
%     %     cp = cvpartition(TimeSeriesGroup,'k',10); % 10-fold stratified cross-validation
%     %     cvMCR = crossval('mcr',TS_DataMat(:,ifeat(i)),TimeSeriesGroup,'predfun',classf,'partition',cp);
% 
%         title(['[' num2str(m_ids(ifeat(i))) ']' mlab{ifeat(i)} ' [' mkw{ifeat(i)} '] -- ' ClassMethod ' = ' num2str(teststat(i))],'interpreter','none')
%         set(gca,'YTick',[])
%     %     if i < i+ ntoplot
%     %         set(gca,'XTick',[]);
%     %     end
%         if i == ntoplot + 1
%             leglabs = cell(Ng*2,1);
%             for j = 1:Ng
%                 leglabs{2*(j-1)+1} = [kwgs{j} ' (' num2str(length(gi{j})) ')'];
%                 leglabs{2*(j-1)+2} = '';
%             end
%             legend(leglabs)
%         end
%     end

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