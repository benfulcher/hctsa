% ------------------------------------------------------------------------------
% Run_IndividualFeatures
% ------------------------------------------------------------------------------
% 
% Script that uses TSQ_IndividualFeatures to search for individual features that
% help distinguish a known classification of the time series, and compares to a
% shuffled distribution.
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

% ------------------------------------------------------------------------------
%% Set Parameters
% ------------------------------------------------------------------------------

% ClassMethod: the classification method to use. The following are possible options:
%            'linclasscv', 'linclass', 'linclasscv', 'knn', 'knncv'
ClassMethod = 'linclass';

% cvalid: in case you choose a cross-validation method
cvalid_0 = {'kfold',10,10}; % 10 repeats! (to reduce variance)
cvalid_shuffle = {'kfold',10,1}; % shuffled should always use just 1 repeat (loops manually)

% NumRepeats: how many times to repeat the calculation using permuted labels
NumRepeats = 5;

% norcl: what local file to load data from
norcl = 'norm';

% Whether to save a record of selected features to file
SaveToFile = 0;


% ------------------------------------------------------------------------------
% Load information from file
% ------------------------------------------------------------------------------
load('HCTSA_N.mat','TimeSeries','Operations');

% ------------------------------------------------------------------------------
% First run TSQ_IndividualFeatures
% ------------------------------------------------------------------------------
nplots = 0;
clear plotopts
plotopts.histogram = 0; % don't plot histogram, we do that below.
plotopts.nfigs = nplots;

[ifeat,testStat,testspread] = TSQ_IndividualFeatures(norcl,ClassMethod,cvalid_0,0,plotopts);

% ------------------------------------------------------------------------------
% Plot histogram of classification errors
% ------------------------------------------------------------------------------
nbins = 20;
figure('color','w');
[N,x] = hist(testStat,nbins);
bar(x,N,'FaceColor','w','EdgeColor','k')
xlabel('Individual classifiction error (%)')
ylabel('Frequency')
titletext = sprintf(['Range = [%4.2f, %4.2f], Mean = %4.2f, Std = %4.2f'], ...
            min(testStat(:)),max(testStat(:)),mean(testStat(~isnan(testStat))),std(testStat(~isnan(testStat))));
title(titletext)
set(gcf,'Position',[400,600,500,270]);

% ------------------------------------------------------------------------------
% List the top operations
% ------------------------------------------------------------------------------
SayWhat = 100;
fn = sprintf('topops_%s',datestr(now,1));
if SaveToFile, fid = fopen([fn '.txt'],'w','n');
else fid = 1; % write to file
end
for i = 1:SayWhat
    if isempty(testspread) % no cross-validation, so no test spread available
        fprintf(fid,'[%u] %s (%s): %4.2f%%\n', ...
                        Operations(ifeat(i)).ID,Operations(ifeat(i)).Name, ...
                        Operations(ifeat(i)).Keywords,testStat(i));
    else
        fprintf(fid,'[%u] %s (%s): %4.2f% +/- %4.2f%%\n', ...
                        Operations(ifeat(i)).ID,Operations(ifeat(i)).Name, ...
                        Operations(ifeat(i)).Keywords,testStat(i),testspread(i));
    end
end

if SaveToFile
    fclose(fid);
    fprintf(1,'Saved %s.txt!!\n',fn)
    
    % Also results as a .mat file
    save([fn '.mat'],'Operations','ifeat','testStat','testspread');
    fprintf(1,'SAVED %s.mat!!!\n',fn)
end

% ------------------------------------------------------------------------------
%% Get randomized individual errors to check robustness
% ------------------------------------------------------------------------------
testStatR = zeros(length(Operations),NumRepeats);
for i = 1:NumRepeats
    tic
    fprintf(1,'REPEAT %u / %u...\n',i,NumRepeats)
    [~,errs] = TSQ_IndividualFeatures(norcl,ClassMethod,cvalid_shuffle,1,plotopts);
    testStatR(:,i) = errs; % just store mean errors each time
    fprintf(1,'REPEAT %u / %u took %s...\n\n',i,NumRepeats,BF_thetime(toc))
end

% ------------------------------------------------------------------------------
% Plot the histogram
% ------------------------------------------------------------------------------
nbins = 20;
figure('color','w');
[N,x] = hist(testStatR(:),nbins);
bar(x,N,'FaceColor','w','EdgeColor','k')
xlabel('Individual Classifiction Error (%)')
ylabel('Frequency')
TitleText = sprintf(['Randomized labeling for %u repeats. Range = [' ...
            '%4.2f, %4.2f], Mean = %4.2f, Std = %4.2f'], ...
            NumRepeats, min(testStatR(:)), max(testStatR(:)), ...
            mean(testStatR(:)), std(testStatR(:)))
title(TitleText)
set(gcf,'Position',[400,600,570,270]);

% ------------------------------------------------------------------------------
% Plot a histogram with both:
% ------------------------------------------------------------------------------
nbins = 20;
figure('color','w');
colormap(BF_getcmap('set1',2,0))
x = linspace(0,max(testStatR(:)),nbins+1);
x = x(1:end-1)+(x(2)-x(1))/2; % centre bins
N1 = hist(testStat(~isnan(testStat)),x)/sum(~isnan(testStat));
N2 = hist(testStatR(~isnan(testStatR)),x)/sum(~isnan(testStatR(:)));
bar(x,[N1;N2]');%,'FaceColor','w','EdgeColor','k');
xlabel(sprintf('%s Classifiction Error (%)',ClassMethod))
ylabel('Frequency')
set(gcf,'Position',[400,600,570,270]);
legend('unshuffled','shuffled')