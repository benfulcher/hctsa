% ------------------------------------------------------------------------------
% TSQ_plot_distribution
% ------------------------------------------------------------------------------
% 
% Input a vector (FeatureVector) and a labeling of its elements (GroupLabels),
% and this will plot a distribution of it.
% 
% Couples nicely with the output of TSQ_IndividualFeatures: e.g., can loop over
% subplots to plot the distributions of the top-performing operations.
% 
%---HISTORY:
% Ben Fulcher, 2014-05-09
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2014,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

function out = TSQ_plot_distribution(FeatureVector,GroupLabels,GroupNames, ...
                                        TitleText,ks_or_hist,MakeNewFigure)

% ------------------------------------------------------------------------------
% Check Inputs:
% ------------------------------------------------------------------------------
if nargin < 3
    GroupNames = {};
end
if nargin < 4
    TitleText = '';
end
if nargin < 5 || isempty(ks_or_hist)
    ks_or_hist = 'ks'; % kernel-smoothed distribution by default (HIGHLY RECOMMENDED!!)
end
if nargin < 6 || isempty(MakeNewFigure)
    MakeNewFigure = 1;
end

% ------------------------------------------------------------------------------
% Preliminaries
% ------------------------------------------------------------------------------

TheLabels = unique(GroupLabels);
NumGroups = length(TheLabels);

My_ColorMap = 'set1';
MarkerColors = BF_getcmap(My_ColorMap,max(3,NumGroups),1);


% ------------------------------------------------------------------------------
% Make the Plot
% ------------------------------------------------------------------------------

if MakeNewFigure
    figure('color','w'); box('on');
end

switch ks_or_hist
case 'ks'
    hold on
    for i = 1:NumGroups
        TheFeatureVectorClassi = FeatureVector(GroupLabels==TheLabels(i));
        [f, x] = ksdensity(TheFeatureVectorClassi);
        plot(x,f,'color',MarkerColors{i},'LineWidth',2);
        % Add dots:
        r = (arrayfun(@(m)find(x >= m,1,'first'),TheFeatureVectorClassi));
        plot(x(r),f(r),'o','MarkerFaceColor',MarkerColors{i},'MarkerEdgeColor',MarkerColors{i})
    end
    
case 'hist'
    % Define bins:
    nbins = 10;
    x = linspace(min(FeatureVector),max(FeatureVector),nbins+1);
    x = x(1:end-1)+(x(2)-x(1))/2; % centre bins
    histwidth = x(2)-x(1);
    
    if histwidth == 0
        error('Feature vector has zero range');
    end

    % Make a histogram for each group:
    Histograms = cell(NumGroups,1);
    for i = 1:NumGroups
        TheFeatureVectorClassi = FeatureVector(GroupLabels==TheLabels(i));
        Histograms{i} = hist(TheFeatureVectorClassi(~isnan(TheFeatureVectorClassi)),x);
    end
    
    % This actually just looks quite ridiculous:
    plotbar = @(z,N,c) rectangle('Position',[z,0,histwidth/2,N],'FaceColor',c);
    for k = 1:length(x);
        for j = 1:NumGroups
            if Histograms{j}(k) > 0
                plotbar(x(k)-histwidth/2,Histograms{j}(k),MarkerColors{j});
            end
        end
    end
    colormap(vertcat(MarkerColors{:}));
end

% ------------------------------------------------------------------------------
% Add a legend
% ------------------------------------------------------------------------------
if ~isempty(GroupNames)
    LegendText = {};
    for i = 1:length(GroupNames),
        LegendText{end+1} = GroupNames{i};
        LegendText{end+1} = '';
    end
    legend(LegendText)
end

ylabel('Probability density')

% ------------------------------------------------------------------------------
% Set the title
% ------------------------------------------------------------------------------
if ischar(TitleText) && ~isempty(TitleText)
    title(TitleText,'interp','none')
end

end
